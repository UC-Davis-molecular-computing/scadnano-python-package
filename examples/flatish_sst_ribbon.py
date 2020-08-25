import argparse

import scadnano as sc


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Create a scadnano design for a flat(ish) SST zig-zag ribbon of a specified width.')
    upper_limit_width = 50
    parser.add_argument('-w', '--width', required=True, type=int,
                        help='width of ribbon (number of tiles in one zig or one zag); '
                             f'must be even number between 4 and {upper_limit_width}',
                        choices=range(4, upper_limit_width + 1, 2))
    parser.add_argument('-d', '--duple', action='store_true',
                        help='whether to use duples (double-tiles) to turn direction; '
                             'otherwise 1.5-tiles are used')
    args = parser.parse_args()
    width: int = args.width
    duple: bool = args.duple
    design = create_design(width=width, duple=duple)
    print(design.to_json())
    design.write_scadnano_file(directory='output_designs')


def create_design(width: int, duple: bool) -> sc.Design:
    helices = [sc.Helix() for _ in range(width + 4)]
    design = sc.Design(helices=helices, grid=sc.Grid.none)

    add_tick_marks(design)
    add_zig_tiles(design, width)
    add_zag_tiles(design, width)
    add_zig_seed_tiles(design, width)
    add_duples(design, width, duple)
    add_seed(design, width)

    return design


def add_tick_marks(design: sc.Design) -> None:
    num_helices = len(design.helices)
    for idx, helix in design.helices.items():
        tick_start = num_helices - idx - 2
        if tick_start % 2 == 1:
            tick_start += 1
        helix.major_tick_start = tick_start
        helix.major_tick_periodic_distances = [12, 9] if idx % 2 == 0 else [10, 11]


def add_zig_tiles(design: sc.Design, width: int) -> None:
    green = sc.Color(hex_string='#009900')
    low_helix_idx = 4
    high_helix_idx = low_helix_idx + width - 3
    offset_5p_even = design.helices[low_helix_idx].major_tick_start + 12 + (width // 2 + 1) * 21
    offset_5p_odd = offset_5p_even - 12
    for idx in range(low_helix_idx, high_helix_idx + 1, 2):
        # 12-9-11-10 strand ("even" strand)
        design.strand(idx, offset_5p_even).move(-12).move(-9) \
            .cross(idx - 1).move(11).move(10).with_color(green)
        offset_5p_even -= 23

        # 11-10-12-9 strand ("odd" strand)
        design.strand(idx + 1, offset_5p_odd).move(-11).move(-10) \
            .cross(idx).move(12).move(9).with_color(green)
        offset_5p_odd -= 23


def add_zag_tiles(design: sc.Design, width: int) -> None:
    red = sc.Color(hex_string='#ff0000')
    low_helix_idx = 3
    high_helix_idx = low_helix_idx + width - 3
    offset_5p_odd = design.helices[low_helix_idx].major_tick_start + (width // 2 + 1) * 21
    offset_5p_even = offset_5p_odd - 11
    for idx in range(low_helix_idx, high_helix_idx + 1, 2):
        # 11-10-12-9 strand ("odd" strand)
        design.strand(idx, offset_5p_odd).move(-11).move(-10) \
            .cross(idx - 1).move(12).move(9).with_color(red)
        offset_5p_odd -= 23

        # 12-9-11-10 strand ("even" strand)
        design.strand(idx + 1, offset_5p_even).move(-12).move(-9) \
            .cross(idx).move(11).move(10).with_color(red)
        offset_5p_even -= 23


def add_zig_seed_tiles(design: sc.Design, width: int) -> None:
    dark_green = sc.Color(hex_string='#006600')
    low_helix_idx = 1
    high_helix_idx = low_helix_idx + width - 2
    offset_5p_odd = design.helices[low_helix_idx].major_tick_start + (width // 2 + 1) * 21
    offset_5p_even = offset_5p_odd - 11

    for idx in range(low_helix_idx, high_helix_idx + 1, 2):
        # 11-10-12-9 strand ("odd" strand)
        design.strand(idx, offset_5p_odd).move(-11).move(-10) \
            .cross(idx - 1).move(12).move(9).with_color(dark_green)
        offset_5p_odd -= 23

        if idx < high_helix_idx:
            # 12-9-11-10 strand ("even" strand)
            design.strand(idx + 1, offset_5p_even).move(-12).move(-9) \
                .cross(idx).move(11).move(10).with_color(dark_green)
            offset_5p_even -= 23


def duple_at(design: sc.Design, color: sc.Color, duple: bool, start_offset: int, top_helix: int) -> None:
    mid_helix = top_helix - 1
    bot_helix = mid_helix - 1
    if duple:
        design.strand(top_helix, start_offset) \
            .move(-11).move(-10) \
            .cross(mid_helix) \
            .move(4).loopout(mid_helix, 4).move(-4) \
            .move(-9) \
            .cross(bot_helix) \
            .move(11).move(10) \
            .cross(mid_helix) \
            .move(-4).loopout(mid_helix, 4).move(4) \
            .move(9) \
            .with_color(color)
    else:
        raise NotImplementedError()


def add_seed_duple(design: sc.Design, width: int, duple: bool, color: sc.Color) -> None:
    top_helix = width + 1
    start_offset = 2 + 21 * 2
    duple_at(design, color, duple, start_offset, top_helix)


def add_zig_to_zag_duple(design: sc.Design, width: int, duple: bool, color: sc.Color) -> None:
    top_helix = width + 3
    start_offset = 21 * 3
    duple_at(design, color, duple, start_offset, top_helix)


def add_zag_to_zig_duple(design: sc.Design, width: int, duple: bool, color: sc.Color) -> None:
    top_helix = 3
    start_offset = design.helices[top_helix].major_tick_start + (width // 2 + 2) * 21
    duple_at(design, color, duple, start_offset, top_helix)


def add_duples(design: sc.Design, width: int, duple: bool) -> None:
    black = sc.Color(hex_string='#000000')
    add_seed_duple(design, width, duple, black)
    add_zig_to_zag_duple(design, width, duple, black)
    add_zag_to_zig_duple(design, width, duple, black)


def add_seed(design: sc.Design, width: int) -> None:
    blue = sc.Color(hex_string='#0066cc')
    start_helix = 0
    start_offset = design.helices[start_helix].major_tick_start + (width // 2 + 1) * 21
    strand_builder = design.strand(start_helix, start_offset).move(-9)
    for helix in range(1, width + 1, 2):
        strand_builder.move(-12).loopout(helix, 8).move(-11)
        if helix < width:
            strand_builder.loopout(helix + 1, 8)
    strand_builder.with_color(blue)


if __name__ == '__main__':
    main()
