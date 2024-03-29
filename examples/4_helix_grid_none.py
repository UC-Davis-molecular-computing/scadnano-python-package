import scadnano as sc


def create_design() -> sc.Design:
    length = 10
    helices = [
        sc.Helix(max_offset=length, position=sc.Position3D(x=0, y=0, z=2.5), roll=0),
        sc.Helix(max_offset=length, position=sc.Position3D(x=3, y=3, z=0), roll=0),
        sc.Helix(max_offset=length, position=sc.Position3D(x=8, y=-3, z=0), roll=0),
        sc.Helix(max_offset=length, position=sc.Position3D(x=11, y=1, z=0), roll=0),
    ]
    stap_ss = sc.Domain(0, True, 0, length)
    scaf_ss = sc.Domain(0, False, 0, length)
    stap = sc.Strand([stap_ss])
    scaf = sc.Strand([scaf_ss], color=sc.default_scaffold_color)
    strands = [stap, scaf]
    design = sc.Design(helices=helices, strands=strands, grid=sc.Grid.none)

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
