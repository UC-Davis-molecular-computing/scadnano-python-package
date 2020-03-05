from dataclasses import dataclass
import scadnano as sc

@dataclass
class Biotin(sc.Modification):
    def __init__(self, location: sc.ModLocation, offset: int = None):
        display = 'B'
        if location == sc.ModLocation.prime5:
            idt = r'/5Biosg/'
            super().__init__(location=location, idt=idt, json_short=idt, display=display)
        elif location == sc.ModLocation.prime3:
            idt = r'/3Bio/'
            super().__init__(location=location, idt=idt, json_short=idt, display=display)
        elif location == sc.ModLocation.internal:
            idt = r'/iBiodT/'
            super().__init__(location=location, offset=offset, idt=idt, json_short=idt, display=display,
                             attached_to_base=True, allowed_bases=['T'])

@dataclass
class Cy3(sc.Modification):
    def __init__(self, location: sc.ModLocation, offset: int = None):
        display = 'Cy3'
        if location == sc.ModLocation.prime5:
            idt = r'/5Cy3/'
            super().__init__(location=location, idt=idt, json_short=idt, display=display)
        elif location == sc.ModLocation.prime3:
            idt = r'/3Cy3Sp/'
            super().__init__(location=location, idt=idt, json_short=idt, display=display)
        elif location == sc.ModLocation.internal:
            idt = r'/iCy3/'
            super().__init__(location=location, offset=offset, idt=idt, json_short=idt, display=display)

@dataclass
class Cy5(sc.Modification):
    def __init__(self, location: sc.ModLocation, offset: int = None):
        display = 'Cy5'
        if location == sc.ModLocation.prime5:
            idt = r'/5Cy5/'
            super().__init__(location=location, idt=idt, json_short=idt, display=display)
        elif location == sc.ModLocation.prime3:
            idt = r'/3Cy5Sp/'
            super().__init__(location=location, idt=idt, json_short=idt, display=display)
        elif location == sc.ModLocation.internal:
            idt = r'/iCy5/'
            super().__init__(location=location, offset=offset, idt=idt, json_short=idt, display=display)
