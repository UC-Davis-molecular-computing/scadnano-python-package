import scadnano as sc
import modifications as mod
import dataclasses

def create_design() -> sc.Design:
    helices = [sc.Helix(max_offset=24) for _ in range(2)]
    design = sc.Design(helices=helices, grid=sc.square)

    design.strand(0,0).move(8).cross(1).move(-8).as_circular()
    design.strand(0,8).move(8).loopout(1, 5).move(-8).as_circular()
    design.strand(0,16).move(8).cross(1).move(-8)

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
