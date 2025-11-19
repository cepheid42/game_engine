from dataclasses import dataclass, field

@dataclass
class EMParams:
    save_interval: int = 1
    nhalo: int = 0
    pml_depth: int = 10
    pml_grade: float = 3.5
    pml_alpha_max: float = 0.2
    em_bcs: tuple = (2, 2, 2, 2, 2, 2)
    
    def __repr__(self):
        return (f'inline constexpr auto em_save_interval = {self.save_interval}zu;\n'
                '\n'
                f'inline constexpr auto PMLDepth    = {self.pml_depth}zu;\n'
                f'inline constexpr auto PMLGrade    = {self.pml_grade};\n'
                f'inline constexpr auto PMLAlphaMax = {self.pml_alpha_max};\n'
                '//inline constexpr auto PMLKappaMax = 1.0;\n'
                '\n'
                f'inline constexpr auto nHalo = {self.nhalo}zu;\n'
                '\n'
                '// Periodic = 0, PML = 1, Reflecting = 2\n'
                f'inline constexpr std::array BCSelect = {{{'zu, '.join(str(i) for i in self.em_bcs) + 'zu'}}};\n'
        )

em = EMParams()
print(em)