from pathlib import Path


def update_header(params, project_path):
    print('Updating header...', end=' ', flush=True)
    print()
    def fancy_header(section_name):
        begin = '/*' + 64 * '-' + '/\n'
        close = '/' + 64 * '-' + '*/\n'
        return begin + f'/-{section_name: ^63}-/\n' + close

    # project_path = Path(project_path)
    # header_path = project_path / "params/"

    open_guard = ('#ifndef PROGRAM_PARAM_HPP\n'
                  '#define PROGRAM_PARAM_HPP\n')
    close_guard = '#endif // PROGRAM_PARAM_HPP'
    includes = ('#include <array>\n'
                '#include <string>\n'
                '#include <tuple>\n')
    print(open_guard + includes + '\n' + fancy_header('Simulation') + close_guard)

update_header(1, 1)
    # program_params = (
    #
    #     '/*---------------------------------------------------------------/\n'
    #     '/-                        EM Parameters                         -/\n'
    #     '/---------------------------------------------------------------*/\n'
    #     'enum class EMFace { X, Y, Z };\n'
    #     'enum class EMSide { Lo, Hi };\n'
    #     '\n'
    #     f'inline constexpr auto em_save_interval = {em_params.save_interval}zu;\n'
    #     '\n'
    #     f'inline constexpr auto PMLDepth    = {em_params.pml_depth}zu;\n'
    #     f'inline constexpr auto PMLGrade    = {em_params.pml_grade};\n'
    #     f'inline constexpr auto PMLAlphaMax = {em_params.pml_alpha_max};\n'
    #     '//inline constexpr auto PMLKappaMax = 1.0;\n'
    #     '\n'
    #     f'inline constexpr auto nHalo = {em_params.nhalo}zu;\n'
    #     '\n'
    #     '// Periodic = 0, PML = 1, Reflecting = 2\n'
    #     f'inline constexpr std::array BCSelect = {{{bc_str}}};\n'
    #     '\n'
    #     '/*---------------------------------------------------------------/\n'
    #     '/-                     Particle Parameters                      -/\n'
    #     '/---------------------------------------------------------------*/\n'
    #     'using collision_spec = std::tuple<std::string, std::string, double, double, int, bool>;\n'
    #     'enum class ParticleBCType { Static, Reflecting, Periodic, Outflow };\n'
    #     '\n'
    #     f'inline constexpr auto particle_save_interval = {particles.save_interval}zu;\n'
    #     f'inline constexpr auto interpolation_order = {particles.interp_order}zu;\n'
    #     '\n'
    #     f'inline constexpr auto PBCSelect = {particle_bcs[particles.particle_bcs]};\n'
    #     '\n'
    #     f'inline constexpr std::array particle_data = {{{particle_data}}};\n'
    #     '\n'
    #     f'inline constexpr std::array<collision_spec, {len(collision_str)}> collision_params = {{'
    #     f'{collisions}'
    #     '};\n'
    #     '\n'
    #     '#endif //PROGRAM_PARAM_HPP\n'
    # )
    #
    # with open(header_path, 'w+') as f:
    #     cur_header = f.read()
    #     if cur_header != program_params:
    #         f.write(program_params)
    # print('Done.')