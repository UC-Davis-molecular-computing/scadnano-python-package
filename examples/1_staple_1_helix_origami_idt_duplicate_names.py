import scadnano as sc

def main():
    ss1_r = sc.Domain(0, True, 0, 4)
    ss2_r = sc.Domain(0, True, 4, 8)
    ss3_r = sc.Domain(0, True, 8, 12)
    ss_l = sc.Domain(0, False, 0, 12)

    s1_r = sc.Strand([ss1_r], idt=sc.IDTFields('s1_r'))
    s2_r = sc.Strand([ss2_r], idt=sc.IDTFields('s1_r'))
    s3_r = sc.Strand([ss3_r], idt=sc.IDTFields('s1_r'))
    s_l = sc.Strand([ss_l], idt=sc.IDTFields('s_l'))

    strands = [s1_r, s2_r, s3_r, s_l]

    design = sc.DNADesign(strands=strands, grid=sc.square)

    design.assign_dna(s_l, 'AGTT'*3)

    return design

if not sc.in_browser() and __name__ == '__main__':
    design = main()
    design.write_idt_bulk_input_file(directory='idt')
