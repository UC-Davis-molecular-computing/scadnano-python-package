import scadnano as sc


def create_design() -> sc.Design:
    helices = [sc.Helix(max_offset=100) for _ in range(4)]
    design = sc.Design(helices=helices, grid=sc.square)

    red = sc.Color(255, 0, 0)
    dark_red = sc.Color(150, 0, 0)
    green = sc.Color(0, 255, 0)
    dark_green = sc.Color(0, 150, 0)
    blue = sc.Color(0, 0, 255)
    dark_blue = sc.Color(0, 0, 150)
    black = sc.Color(0, 0, 0)

    design.draw_strand(0, 0) \
        .extension_5p(num_bases=5).with_domain_color(red) \
        .move(8).with_domain_color(green) \
        .loopout(1, 5).with_domain_color(dark_blue) \
        .move(-8).with_domain_color(dark_red) \
        .cross(2) \
        .move(8).with_domain_color(dark_green) \
        .cross(3) \
        .move(-8) \
        .extension_3p(num_bases=5).with_domain_color(black) \
        .with_color(blue)

    return design

if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
