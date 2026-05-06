from dataclasses import dataclass

@dataclass
class EMParams:
    save_interval: int = 10
    nhalo: int = 0
    pml_depth: int = 10
    pml_grade: float = 3.5
    pml_alpha_max: float = 0.2
    em_bcs: tuple = (2, 2, 2, 2, 2, 2)
    applied_fields: str = ''