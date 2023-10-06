from typing import List
import math
import numpy as np
from parameters import *

class InflowPressureRelation():
    def __init__(self, *, j_index, pressure_reservoir, pressure_bubble, watercut = 0):
        self.j = j_index
        self.p_res = pressure_reservoir
        self.p_bubble = pressure_bubble
        self.watercut = watercut
        self.qmax = (self.j * (self.p_res - self.p_bubble) /
                     (1 - 0.2 * (self.p_bubble / self.p_res) - 0.8 * (self.p_bubble / self.p_res) ** 2))
        self.qb = j_index * (pressure_reservoir - pressure_bubble)
        self.qmax = self.qb + j_index * pressure_bubble / 1.8
        self.p_wfg = watercut * (pressure_reservoir - self.qmax / j_index)

    def composite_ipr(self, p):
        p_b = self.p_bubble
        p_r = self.p_res
        p_wfg = self.p_wfg
        qmax = self.qmax
        qb = self.qb
        fw = self.watercut
        fo = 1 - fw
        j = self.j

        if p >= p_b:
            return j * (p_r - p)
        elif (p < p_b) and (p > p_wfg):
            A = (p + 0.125 * fo * p_b - fw * p_r) / (0.125 * fo * p_b)
            B = fw / (0.125 * fo * p_b * j)
            C = 2 * A * B + 80 / (qmax - qb)
            D = A**2 - 80 * qb / (qmax - qb) - 81

            if B == 0:
                return - D / C
            else:
                return ( - C + math.sqrt(C ** 2 - 4 * D * B ** 2)) / (2 * B ** 2)
        else:
            CD = fw * 0.001 * qmax / j + fo * 0.125 * p_b * (math.sqrt(81 - 80 * (0.999 * qmax - qb) / (qmax - qb)) - 1)  
            CG = 0.001 * qmax
            beta = CD / CG
            return (p_wfg + qmax * beta - p) / beta

    def ipr(self, q) -> float:
        p= self.p_res - q / self.j
        if p >= self.p_bubble:
            return p
        else:
            x=(-0.2 + math.sqrt(0.04 - 3.2 * (q / self.qmax - 1))) / 1.6
            return x * self.p_res
        
    def ipr_list(self) -> List[float]:
        #flow_values = range(1, int(self.qmax) + 1)
        #pressure_values = [self.ipr(q) for q in flow_values]
        flow_values = [self.composite_ipr(p) for p in range(1, self.p_res, 5)]
        return flow_values, range(1, self.p_res, 5)

class VerticalLiftPerformance():

    def __init__(self, *, sg_oil, sg_gas, gas_oil_ratio, oil_volume_factor, pressure_reservoir,
                temperature_reservoir, diameter, liquid_viscosity, gas_viscosity, depth, pipe_length=None, 
                thp = 30, tht = 77, watercut = 0, ed = 4.572e-5, sigma_oil = 1, water_viscosity = 0.64, 
                sigma_water = 48.1, sg_water = 1.04, water_volume_factor = 1.01):
        if pipe_length ==None:
            pipe_length = depth
        self.sg_liquid = sg_oil
        self.sg_gas = sg_gas
        self.gor = gas_oil_ratio
        self.b_oil = oil_volume_factor
        self.pressure = pressure_reservoir
        self.temp = temperature_reservoir
        self.tht = tht
        self.thermalgradient = (temperature_reservoir - tht) / depth
        self.diameter = diameter
        self.oil_visco = liquid_viscosity
        self.gas_visco = gas_viscosity
        self.length = pipe_length
        self.ed = ed
        self.sigma_oil = sigma_oil
        self.glr = gas_oil_ratio * (1 - watercut)
        self.water = watercut
        self.water_visco = water_viscosity
        self.sigma_water = sg_water
        self.water_sg = sg_water
        self.b_water = water_volume_factor
        self.depth = depth
        self.thp = thp
        if watercut != 0:
            self.oilcut = 1 - watercut
        else:
            self.oilcut = 1

    def friction_factor(self, Re):
        return 0.25 / (math.log10(self.ed / self.diameter / 3.7 + 5.74 / (Re ** 0.9))) ** 2
                
    def z_factor(self, P, t):
        y = self.sg_gas
        #t = 1.8 * (self.temp - 273) + 492

        p_c = 756.8 - 131 * y - 3.6 * y ** 2 #psia
        t_c = 169.2 + 349.5 * y - 74 * y ** 2 #rankine

        Tr = t / t_c
        Pr = P / p_c

        a = 1.39 * (Tr - 0.92) ** 0.5 - 0.36 * Tr - 0.101
        b = (0.62 - 0.23 * Tr) * Pr + (0.066 / (Tr - 0.86) - 0.037) * Pr ** 2 + 0.32 * Pr ** 6 / (10** (9 * (Tr - 1)))
        c = (0.132 - 0.32 * math.log10(Tr))
        d = 10 ** (0.3106 - 0.49 * Tr + 0.1824 * Tr **2)
        zfact = a + (1 - a) * math.exp(-b) + c * Pr **d

        return zfact

    def hold_up(self, rho_l, liquid_visco, sigma_l, q, z, pressure, temp):
    
        A = 3.1415 * (self.diameter) ** 2 / 4
        Nl = 0.15726 * liquid_visco * (rho_l * sigma_l ** 3) ** (-0.25)
        CNl = 0.061 * Nl ** 3 - 0.0929 * Nl ** 2 + 0.0505 * Nl + 0.0019
        
        vs_l = 5.616 * q / (86400 * A) * (self.b_oil * (1 - self.water) + self.b_water * self.water)
        vs_g = q * self.glr / 86400 / A * 14.7 / pressure * temp / 520 * z
        
        Nlv = 1.938 * vs_l * (rho_l / sigma_l) ** (0.25)
        Ngv = 1.938 * vs_g * (rho_l / sigma_l) ** (0.25)
        Nd = 120.872 * self.diameter * (rho_l / sigma_l) ** (0.5)

        H = Nlv / (Ngv ** 0.575) * (pressure / 14.7) ** 0.1 * CNl / Nd 
        
        Hlv = math.sqrt((0.0047 + 1123.32 * H + 729489.64 * H * H) / (1 + 1097.1566 * H + 722153.97 * H * H))
        B = Ngv * Nlv ** 0.38 / (Nd ** 2.14)

        if B <= 0.025:
            fi = 27170 * B ** 3 - 317.52 * B ** 2 + 0.5472 * B + 0.9999
        elif B > 0.055:
            fi = 2.5714 * B + 1.5962
        else:
            fi = -533.33 * B ** 2 + 58.524 * B + 0.1171

        return Hlv * fi, vs_g, vs_l
         

    def pressure_gradient(self, q, pressure, temp):
        temp = temp + 492
        M = 350.52 * self.sg_liquid * (1 - self.water)  + self.sg_gas * 0.0764 * self.gor /(1 - self.water)
        
        z = self.z_factor(pressure, temp)
        rho_l = (62.4 * self.sg_liquid + self.gor * 0.0764 * self.sg_gas / 5.614)  * (1 - self.water) + 62.4 * self.water_sg * self.water
        rho_g = 28.967 * self.sg_gas * pressure / z / temp / 10.732

        liquid_visco = self.oil_visco * (1-self.water) + self.water_visco * self.water
        sigma_l = self.sigma_oil * (1-self.water) + self.sigma_water * self.water

        x = self.hold_up(rho_l, liquid_visco, sigma_l, q, z, self.pressure, temp)
        h = x[0]
        #x = self.hold_up_2(q, rho_l, liquid_visco, sigma_l, pressure, temp, z)
        #h = x[0]

        Re = 0.022 * q * M / self.diameter / (liquid_visco ** h) / (self.gas_visco ** (1 - h))

        rho_m = rho_l * h + rho_g * (1 - h)
        
        if Re > 2000:
            fric = self.friction_factor(Re)
        else:
            fric = 64 / Re
        
        hydro = rho_m
        #friction = fric * (q * M) ** 2 / (2.9652 * 10 ** 11 * (self.diameter ** 5) * rho_l)
        
        friction = 2 * fric * rho_l * (x[2] / x[0]) ** 2 / (self.diameter / 12) / 144

        Ek = 0 #(x[1] + x[2]) * x[1] * rho_m / 32.17 / pressure / 144
        return (hydro+ friction) / 144 / (1 - Ek)
    
    def pressure_traverse(self, flow_rate):
        depths = np.linspace(0, self.depth, 50)
        temp = self.tht
        p=[]
        dpdz=[]
    
        for i in range(50):

            if i==0:
                p.append(self.thp)
            else:
                dz = (depths[i]-depths[i-1])
                pressure = p[i-1] + dz * dpdz[i-1]
                temp = temp + self.thermalgradient * dz
                p.append(pressure)

            dpdz_step = self.pressure_gradient(flow_rate, p[i], temp)
             
            dpdz.append(dpdz_step)
            
        return p, dpdz

    def vlp_values(self, irp):
        pressure_values = []

        for _ in irp[0]:
            pressure = self.pressure_traverse(_)
            pressure_values.append(pressure[0][-1])

        return irp[0], pressure_values
