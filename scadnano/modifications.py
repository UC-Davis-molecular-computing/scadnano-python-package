from scadnano import Modification5Prime, Modification3Prime, ModificationInternal

biotin_5p = Modification5Prime(display_text='B', id='/5Biosg/', idt_text='/5Biosg/')
biotin_3p = Modification3Prime(display_text='B', id='/3Bio/', idt_text='/3Bio/')
biotin_int = ModificationInternal(display_text='B', id='/iBiodT/', idt_text='/iBiodT/',
                                  allowed_bases=frozenset('T'))

cy3_5p = Modification5Prime(display_text='Cy3', id='/5Cy3/', idt_text='/5Cy3/')
cy3_3p = Modification3Prime(display_text='Cy3', id='/3Cy3Sp/', idt_text='/3Cy3Sp/')
cy3_int = ModificationInternal(display_text='Cy3', id='/iCy3/', idt_text='/iCy3/')

cy5_5p = Modification5Prime(display_text='Cy5', id='/5Cy5/', idt_text='/5Cy5/')
cy5_3p = Modification3Prime(display_text='Cy5', id='/3Cy5Sp/', idt_text='/3Cy5Sp/')
cy5_int = ModificationInternal(display_text='Cy5', id='/iCy5/', idt_text='/iCy5/')
