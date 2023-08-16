import scadnano as sc
import modifications as mod
import dataclasses


def create_design() -> sc.Design:
    # '''
    #   0123456789012345678901234567890123456789
    # 0 [---+[--------+[----------+
    #       |         |           |
    # 1 [---+<--------+<----------+
    #
    # angle (fraction of 360)
    #      5/10.5
    #              (15-10.5)/10.5 = 4.5/10.5
    #                        (27-21)/10.5 = 6/10.5
    # '''
    # design2h = sc.Design(helices=[sc.Helix(max_offset=50) for _ in range(2)], grid=sc.square)
    # # helix 0 forward
    # design2h.draw_strand(0, 0).move(5).cross(1).move(-5)
    # design2h.draw_strand(0, 5).move(10).cross(1).move(-10)
    # design2h.draw_strand(0, 15).move(12).cross(1).move(-12)
    #
    # for helix in design2h.helices.values():
    #     helix.major_ticks = [0, 5, 15, 27]
    #
    # design2h.relax_helix_rolls()

    '''
      0         1         2         3         4         5         6
      012345678901234567890123456789012345678901234567890123456789
    0 [---+[--------+[----------+[------+[--------+[--------+
          |         |           |       |         |         |
    1 [---+<--------+<----------+       |         |         |
                                        |         |         |
    2                            <------+<--------+<--------+

    angle (fraction of 360)
         5/10.5
                 (15-10.5)/10.5 = 4.5/10.5
                           (27-21)/10.5 = 6/10.5
    '''
    helices = [sc.Helix(max_offset=60) for _ in range(3)]
    helices[2].grid_position = (1, 0)
    design3h = sc.Design(helices=helices, grid=sc.square)

    # helix 0 forward
    design3h.draw_strand(0, 0).move(5).cross(1).move(-5)
    design3h.draw_strand(0, 5).move(10).cross(1).move(-10)
    design3h.draw_strand(0, 15).move(12).cross(1).move(-12)
    design3h.draw_strand(0, 27).move(7).cross(2).move(-7)
    design3h.draw_strand(0, 34).move(10).cross(2).move(-10)
    design3h.draw_strand(0, 44).move(10).cross(2).move(-10)

    for helix in design3h.helices.values():
        helix.major_ticks = [0, 5, 15, 27, 34, 44, 54]

    design3h.relax_helix_rolls()

    # return design3h

    helices = [sc.Helix(max_offset=60) for _ in range(3)]
    helices[2].grid_position = (1, 0)
    design3h2 = sc.Design(helices=helices, grid=sc.square)
    design3h2.draw_strand(0, 0).move(5).cross(1).move(-5)
    design3h2.draw_strand(0, 5).move(8).cross(2).move(-8)

    for helix in design3h2.helices.values():
        helix.major_ticks = [0, 5, 13]

    design3h2.relax_helix_rolls()

    return design3h2


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
