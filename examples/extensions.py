import scadnano as sc


def create_design() -> sc.Design:
    helices = [sc.Helix(max_offset=16) for _ in range(3)]
    design = sc.Design(helices=helices, grid=sc.square)

    design.draw_strand(0, 0).extension_5p(5, display_length=2, display_angle=30)\
        .move(16).cross(1).move(-16).loopout(2, 3).move(16).extension_3p(7).with_domain_name("ext_3p")

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
