from scadnano import Modification5Prime, Modification3Prime, ModificationInternal

biotin_5p = Modification5Prime(display_text='B', idt_text='/5Biosg/')
biotin_3p = Modification3Prime(display_text='B', idt_text='/3Bio/')
biotin_int = ModificationInternal(display_text='B', idt_text='/iBiodT/', allowed_bases=frozenset('T'))

cy3_5p = Modification5Prime(display_text='Cy3', idt_text='/5Cy3/')
cy3_3p = Modification3Prime(display_text='Cy3', idt_text='/3Cy3Sp/')
cy3_int = ModificationInternal(display_text='Cy3', idt_text='/iCy3/')

cy5_5p = Modification5Prime(display_text='Cy5', idt_text='/5Cy5/')
cy5_3p = Modification3Prime(display_text='Cy5', idt_text='/3Cy5Sp/')
cy5_int = ModificationInternal(display_text='Cy5', idt_text='/iCy5/')

fam_5p = Modification5Prime(display_text='F', idt_text='/56-FAM/')
fam_3p = Modification3Prime(display_text='F', idt_text='/36-FAM/')

rox_5p = Modification5Prime(display_text='R', idt_text='/56-ROXN/')
rox_3p = Modification3Prime(display_text='R', idt_text='/3Rox_N/')

fluorescein_5p = Modification5Prime(display_text='F', idt_text='/5FluorT/')
fluorescein_3p = Modification3Prime(display_text='F', idt_text='/3FluorT/')
fluorescein_int = ModificationInternal(display_text='F', idt_text='/iFluorT/')