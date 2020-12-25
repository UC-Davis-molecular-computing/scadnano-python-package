import scadnano as sc


def create_design() -> sc.Design:
    helices = []
    num_helices, max_bases = 30, 2000
    for _ in range(num_helices):
        helices.append(sc.Helix(max_bases))
    design = sc.Design(helices=helices, strands=[], grid=sc.square)

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
