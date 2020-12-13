import scadnano as sc
import modifications as mod
import dataclasses

def create_design():
    helices = [sc.Helix(max_offset=30) for _ in range(2)]
    design = sc.Design(helices=helices, grid=sc.square)

    design.strand(0,0).move(10).cross(1).move(-10).as_circular()
    design.strand(0,10).move(10).loopout(1, 5).move(-10).as_circular()
    design.strand(0,20).move(10).cross(1).move(-10)

    return design


if __name__ == '__main__':
    design = create_design()
    design.write_scadnano_file(directory='output_designs')
