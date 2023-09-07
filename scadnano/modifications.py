from scadnano import Modification5Prime, Modification3Prime, ModificationInternal

biotin_5p = Modification5Prime(display_text='B', vendor_code='/5Biosg/')
biotin_3p = Modification3Prime(display_text='B', vendor_code='/3Bio/')
biotin_int = ModificationInternal(display_text='B', vendor_code='/iBiodT/', allowed_bases=frozenset('T'))

cy3_5p = Modification5Prime(display_text='Cy3', vendor_code='/5Cy3/')
cy3_3p = Modification3Prime(display_text='Cy3', vendor_code='/3Cy3Sp/')
cy3_int = ModificationInternal(display_text='Cy3', vendor_code='/iCy3/')

cy5_5p = Modification5Prime(display_text='Cy5', vendor_code='/5Cy5/')
cy5_3p = Modification3Prime(display_text='Cy5', vendor_code='/3Cy5Sp/')
cy5_int = ModificationInternal(display_text='Cy5', vendor_code='/iCy5/')

fam_5p = Modification5Prime(display_text='F', vendor_code='/56-FAM/')
fam_3p = Modification3Prime(display_text='F', vendor_code='/36-FAM/')

rox_5p = Modification5Prime(display_text='R', vendor_code='/56-ROXN/')
rox_3p = Modification3Prime(display_text='R', vendor_code='/3Rox_N/')

fluorescein_5p = Modification5Prime(display_text='F', vendor_code='/5FluorT/')
fluorescein_3p = Modification3Prime(display_text='F', vendor_code='/3FluorT/')
fluorescein_int = ModificationInternal(display_text='F', vendor_code='/iFluorT/')
