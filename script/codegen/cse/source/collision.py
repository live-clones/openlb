from sympy import *
from sympy.codegen.ast import CodeBlock, Assignment
from sympy.utilities.iterables import numbered_symbols
from source.cse_utils import SymbolGenerator, CodeBlockPrinter, code, cse
import re

def collision_cse(inputFile):
    # Load expression tree from ouput file
    with open(inputFile, "r") as file:
        load_expressions = file.read()
    results = { }
    exec(load_expressions, globals(), results)
    collisionO = results['collisionO']
    cell_assignments = results['cell_assignments']
    params_symbols = results['params_symbols']
    return_assignments = results['return_assignments']

    # Retrieve dynamics tuple information
    collisionO = collisionO.replace("Expr", "T")
    collisionO = re.sub(r'(descriptors::D\dQ\d{1,2})<.*?>', r'\1<FIELDS...>', collisionO)

    # Combine all assignments
    assignments = cell_assignments
    optional_symbols = params_symbols

    # Current limitations of sympy assignments
    generator = iter(SymbolGenerator('x'))
    cell_backsubstitutions = [ ]
    sub_assignments = [ ]
    for assignment in assignments:
        alias = next(generator)
        cell_backsubstitutions.append((alias, assignment.lhs))
        sub_assignments.append(Assignment(alias, assignment.rhs))

    # Creating aliases and assignments for parameter calls
    optional_assignments = [ ]
    optional_substitutions = [ ]
    for symbol in optional_symbols:
        alias = next(generator)
        optional_assignments.append(Assignment(alias, symbol))
        optional_substitutions.append((symbol, alias))

    # Create code block to perform cse
    block = CodeBlock(*sub_assignments,*return_assignments)
    block = block.subs(optional_substitutions)

    # Perform cse
    block = cse(block, generator)

    # These are required due to current CSE limititations
    substitutions = [ ]
    post_collision = [ ]
    simple_post_collision = True
    for assignment in block.args:
        if assignment.lhs in cell_backsubstitutions[0]:
            atoms = assignment.rhs.free_symbols
            for (lhs, rhs) in cell_backsubstitutions:
                if lhs in atoms or rhs in atoms:
                    simple_post_collision = False
    if simple_post_collision:
        substitutions = cell_backsubstitutions
    else:
        for subs in cell_backsubstitutions:
            post_collision.append(Assignment(subs[1], subs[0]))

    # Apply back-substitution of placeholders
    block = block.subs(substitutions)

    # Detect which optional assignments are required to close the symbolic collision
    optional = [ ]
    for symbol in block.free_symbols:
        if symbol in { assgn.lhs for assgn in optional_assignments }:
            for assgn in optional_assignments:
                if assgn.lhs is symbol:
                    optional.append(assgn)

    # return cse optimized expression
    return '\n'.join([
        "template <typename T, typename... FIELDS>",
        f"""struct CSE<{ collisionO }> {{""",
        "template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>",
        "CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {",
        code(CodeBlock(*optional, *block.args[:-2], *post_collision), symbols=generator.symbols),
        "return { %s, %s };" % (code(block.args[-2].rhs), code(block.args[-1].rhs)),
        "}\n};"
    ])
