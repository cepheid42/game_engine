from dataclasses import dataclass
from enum import Enum

from ..params.header_utils import *


class EMBCType(Enum):
    Periodic = 0
    PML = 1
    Reflecting = 2


@dataclass
class EMParams:
    save_interval: int = 10
    nhalo: int = 0
    pml_depth: int = 10
    pml_grade: float = 3.5
    pml_alpha_max: float = 0.2
    em_bcs: tuple = (2, 2, 2, 2, 2, 2)
    external_fields_file: str = ''
    laser_enabled: bool = False # todo: these should be specs instead
    raman_enabled: bool = False
    rmf_params: tuple = None

    def __repr__(self):
        bc_string = f'{{ {", ".join(self.em_bcs)} }}'

        rmf_enabled = False
        rmf_params_str = ''
        if self.rmf_params is not None:
            rmf_enabled = True
            rmf_params_str = f'{{ {", ".join([str(val) for val in self.rmf_params])} }}'

        return str(
            section_label('EM Parameters') +
            '\n' +
            enum_declaration('EMFace', ['X', 'Y', 'Z']) +
            enum_declaration('EMSide', ['Lo', 'Hi']) +
            '\n' +
            # constexpr_declaration('em_save_interval', f'{int(self.save_interval)}zu') +
            # constexpr_declaration('em_subcycles', f'{int(num_subcycles_em)}zu') +
            # constexpr_declaration('dt_em', float(dt_em)) +
            # '\n' +
            constexpr_declaration('PMLDepth', f'{int(self.pml_depth)}zu') +
            constexpr_declaration('PMLGrade', float(self.pml_grade)) +
            constexpr_declaration('PMLAlphaMax', float(self.pml_alpha_max)) +
            constexpr_declaration('PMLKappaMax', 1.0) +
            '\n' +
            constexpr_declaration('nHalo', f'{int(self.nhalo)}zu') +
            '\n' +
            comment('Periodic = 0, PML = 1, Reflecting = 2') +
            constexpr_declaration('BCSelect', bc_string, typestr='std::array') +
            '\n' +
            constexpr_declaration('laser_enabled', bool2str(self.laser_enabled)) +
            constexpr_declaration('raman_enabled', bool2str(self.raman_enabled)) +
            constexpr_declaration('rmf_enabled', bool2str(rmf_enabled)) +
            '\n' +
            constexpr_declaration('rmf_params', rmf_params_str, typestr='std::array') +
            '\n' +
            constexpr_declaration('external_fields_file', f'"{self.external_fields_file}"sv') +
            '\n'
        )