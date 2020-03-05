from scadnano import Modification5Prime, Modification3Prime, ModificationInternal

biotin_5p = Modification5Prime(display_text='B', json_short='/5Biosg/', idt_text='/5Biosg/')
biotin_3p = Modification3Prime(display_text='B', json_short='/3Bio/', idt_text='/3Bio/')
biotin_int = ModificationInternal(display_text='B', json_short='/iBiodT/', idt_text='/iBiodT/',
                                  allowed_bases=frozenset('T'))

cy3_5p = Modification5Prime(display_text='Cy3', json_short='/5Cy3/', idt_text='/5Cy3/')
cy3_3p = Modification3Prime(display_text='Cy3', json_short='/3Cy3Sp/', idt_text='/3Cy3Sp/')
cy3_int = ModificationInternal(display_text='Cy3', json_short='/iCy3/', idt_text='/iCy3/')

cy5_5p = Modification5Prime(display_text='Cy5', json_short='/5Cy5/', idt_text='/5Cy5/')
cy5_3p = Modification3Prime(display_text='Cy5', json_short='/3Cy5Sp/', idt_text='/3Cy5Sp/')
cy5_int = ModificationInternal(display_text='Cy5', json_short='/iCy5/', idt_text='/iCy5/')

# class Biotin(sc.Modification):
#     def __init__(self, location: sc.ModLocation, offset: int = None):
#         display = 'B'
#         if location == sc.ModLocation.prime5:
#             idt = r'/5Biosg/'
#             super().__init__(location=location, idt=idt, json_short=idt, display=display)
#         elif location == sc.ModLocation.prime3:
#             idt = r'/3Bio/'
#             super().__init__(location=location, idt=idt, json_short=idt, display=display)
#         elif location == sc.ModLocation.internal:
#             idt = r'/iBiodT/'
#             super().__init__(location=location, offset=offset, idt=idt, json_short=idt, display=display,
#                              attached_to_base=True, allowed_bases=['T'])
#
#
# class Cy3(sc.Modification):
#     def __init__(self, location: sc.ModLocation, offset: int = None):
#         display = 'Cy3'
#         if location == sc.ModLocation.prime5:
#             idt = r'/5Cy3/'
#             super().__init__(location=location, idt=idt, json_short=idt, display=display)
#         elif location == sc.ModLocation.prime3:
#             idt = r'/3Cy3Sp/'
#             super().__init__(location=location, idt=idt, json_short=idt, display=display)
#         elif location == sc.ModLocation.internal:
#             idt = r'/iCy3/'
#             super().__init__(location=location, offset=offset, idt=idt, json_short=idt, display=display)
#
#
# class Cy5(sc.Modification):
#     def __init__(self, location: sc.ModLocation, offset: int = None):
#         display = 'Cy5'
#         if location == sc.ModLocation.prime5:
#             idt = r'/5Cy5/'
#             super().__init__(location=location, idt=idt, json_short=idt, display=display)
#         elif location == sc.ModLocation.prime3:
#             idt = r'/3Cy5Sp/'
#             super().__init__(location=location, idt=idt, json_short=idt, display=display)
#         elif location == sc.ModLocation.internal:
#             idt = r'/iCy5/'
#             super().__init__(location=location, offset=offset, idt=idt, json_short=idt, display=display)
