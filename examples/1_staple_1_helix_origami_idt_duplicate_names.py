import scadnano as sc


def create_design() -> sc.Design:
    ss1_r = sc.Domain(0, True, 0, 4)
    ss2_r = sc.Domain(0, True, 4, 8)
    ss3_r = sc.Domain(0, True, 8, 12)
    ss_l = sc.Domain(0, False, 0, 12)

    s1_r = sc.Strand([ss1_r], vendor_fields=sc.VendorFields(), name='s1_r')
    s2_r = sc.Strand([ss2_r], vendor_fields=sc.VendorFields(), name='s1_r')
    s3_r = sc.Strand([ss3_r], vendor_fields=sc.VendorFields(), name='s1_r')
    s_l = sc.Strand([ss_l], vendor_fields=sc.VendorFields(), name='s_l')

    strands = [s1_r, s2_r, s3_r, s_l]

    design = sc.Design(strands=strands, grid=sc.square)

    design.assign_dna(s_l, 'AGTT' * 3)

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_idt_bulk_input_file(directory='idt')
