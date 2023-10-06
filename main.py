import NodalAnalysis
import Visual

def main():
    x, y = NodalAnalysis.curve_analysis()
    Visual.nodal_graph(x, y)

if __name__ == "__main__":
    main()