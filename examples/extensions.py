import scadnano as sc


def main() -> None:
    d = create_design()
    d.write_scadnano_file(directory='output_designs')


def create_design() -> sc.Design:
    width = 8
    helices = [sc.Helix(max_offset=32) for _ in range(3)]
    design = sc.Design(helices=helices, grid=sc.square)

    design.draw_strand(0, 0).extension_5p(5, display_length=2.5, display_angle=45) \
        .with_domain_name("ext_5p 1") \
        .move(width).with_domain_name('domain 1') \
        .cross(1).move(-width) \
        .loopout(2, 3).with_domain_name('loopout 1') \
        .move(width) \
        .extension_3p(7).with_domain_name("ext_3p 1") \
        .with_name('strand1') \
        .with_sequence('T' * 5 + 'G' * width + 'C' * width + 'T' * 3 + 'G' * width + 'A' * 7) \
        .with_modification_5p(sc.cy5_5p) \
        .with_modification_3p(sc.cy3_3p)

    design.draw_strand(0, 24).extension_5p(5, display_length=3.5, display_angle=60) \
        .with_domain_name("ext_5p 2") \
        .move(-width).cross(1).move(width).loopout(2, 3).move(-width) \
        .extension_3p(7).with_domain_name("ext_3p 2") \
        .with_name('strand2') \
        .with_sequence('A' * 5 + 'C' * width + 'G' * width + 'A' * 3 + 'C' * width + 'T' * 7) \

    return design


if __name__ == '__main__':
    main()
