import scadnano as sc
import modifications as mod
import dataclasses


def create_design() -> sc.Design:
    '''
      0         1         2
      012345678901234567890
    0        <--+c+--]
    '''
    helices = [sc.Helix(max_offset=20) for _ in range(1)]
    design = sc.Design(helices=helices, grid=sc.square)
    sb = design.draw_strand(0, 16)
    sb.move(-4)
    sb.cross(0, 11)
    sb.move(-4)

    for helix in design.helices.values():
        helix.major_ticks = [0, 5, 10, 15, 20]

    design.relax_helix_rolls()

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
