import scadnano as sc


def create_design():
    num_domains = 5
    helices = [sc.Helix(max_offset=num_domains * 8)]
    design = sc.Design(helices=helices, strands=[], grid=sc.square)
    design.strand(0, num_domains * 8) \
        .move(-8).with_domain_name('domain 5*') \
        .move(-8).with_domain_name('domain 4*') \
        .move(-8).with_domain_name('domain 3*') \
        .move(-8).with_domain_name('domain 2*') \
        .move(-8).with_domain_name('domain 1*') \
        .with_name('bottom strand')

    # domain names match
    design.strand(0, 0) .move(8).with_domain_name('domain 1') .with_name('top strand 1')

    # domains are aligns but names mismatch
    design.strand(0, 8) .move(8).with_domain_name('domain2') .with_name('top strand 2')

    # domain names match
    design.strand(0, 16) .move(8).with_domain_name('domain 3') .with_name('top strand 3')

    # domain names match but domains are mis-aligned
    design.strand(0, 24) .move(9).with_domain_name('domain 4') .with_name('top strand 4')

    # domain names match but domains are mis-aligned
    design.strand(0, 33) .move(7).with_domain_name('domain 5') .with_name('top strand 5')

    return design


if __name__ == '__main__':
    design = create_design()
    design.write_scadnano_file(directory='output_designs')
