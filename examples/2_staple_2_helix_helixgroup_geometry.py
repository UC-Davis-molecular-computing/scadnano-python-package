import scadnano as sc


def create_design() -> sc.Design:
    group0 = sc.HelixGroup(grid=sc.square)
    group1 = sc.HelixGroup(grid=sc.square, geometry=sc.Geometry(bases_per_turn=18),
                           position=sc.Position3D(0, 3, 0))
    groups = {"group 0": group0, "group 1": group1}
    helices = [sc.Helix(idx=idx, max_offset=40, group=group) for idx, group in
               [(0, "group 0"), (1, "group 1")]]
    design = sc.Design(helices=helices, groups=groups, strands=[])
    design.draw_strand(0, 0).move(40)
    design.draw_strand(0, 40).move(-40)
    design.draw_strand(1, 0).move(40)
    design.draw_strand(1, 40).move(-40)

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
    d.from_scadnano_file('output_designs/2_staple_2_helix_helixgroup_geometry.sc')
    print(f'design: {d.to_json()}')
