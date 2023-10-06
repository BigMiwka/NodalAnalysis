from PressureFlowCurves import VerticalLiftPerformance, InflowPressureRelation
from parameters import *

def intersection_point(irp, vlp):
    result = []

    for i in range(len(irp[1])):
        result.append(abs(irp[1][i] - vlp[1][i]))

    index = result.index(min(result))

    return irp[0][index] ,irp[1][index]

def curve_analysis():
    x = VerticalLiftPerformance(sg_oil = sg_o, sg_gas = sg_g, gas_oil_ratio = gor, oil_volume_factor = B0, pressure_reservoir = p_r, 
                                temperature_reservoir = temp, diameter = d, liquid_viscosity = mu_o, gas_viscosity = mu_g, depth = h)

    y = InflowPressureRelation(j_index=j, pressure_reservoir= p_r, pressure_bubble=p_b)

    irp = y.ipr_list()
    vlp = x.vlp_values(irp)

    return irp, vlp


