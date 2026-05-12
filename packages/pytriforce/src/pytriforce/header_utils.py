class NameGetter:
    def get_name(self):
        return str(self).split(':')[-1]

def constexpr_declaration(name, value, typestr='auto'):
    return f'inline constexpr {typestr} {name} = {value};\n'

def enum_declaration(name, values):
    return f'enum class {name} {{{', '.join(v.get_name() if hasattr(v, 'get_name') else v for v in values)} }};\n'

def section_label(title):
    return '/*' + 46 * '-' + '-/\n/-' + f'{title:^46}' + '-/\n/-' + 46 * '-' + '*/\n'

def comment(msg):
    return f'// {msg}\n'

def tuple_to_array(vals):
    result = f'{{{', '.join(str(v) for v in vals)}}}'
    if len(result) > 80:
        result.replace('{', '{\n\t')
        result.replace(', ', ',\n\t')
        result.replace('}', '\n}')
    return result

def bool2str(val):
    return str(val).lower()
