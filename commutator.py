import sympy as sp

H, A, B, C, D = sp.symbols('H A B C D', commutative=False)

def comm(A, B):
    return A*B - B*A

def simpl(expr):
    return sp.expand(expr)

def ABpow(term, A, B):
    factors = term.as_ordered_factors()
    a_count = 0
    b_count = 0
    for factor in factors:
        if factor.has(A):
            if factor.is_Pow:  
                a_count += factor.exp
            else:
                a_count += 1 
        elif factor.has(B):
            if factor.is_Pow: 
                b_count += factor.exp
            else:
                b_count += 1  
    return a_count, b_count

def filter(expr, H, A, B):
    terms = expr.as_ordered_terms()
    filtered_terms = []
    
    for term in terms:
        factors = term.as_ordered_factors()
        
        for factor in factors:
            if not factor.is_number:  
                if factor == H:
                    a_count, b_count = ABpow(term, A, B)
                    if a_count + b_count <= 4: 
                        filtered_terms.append(term)
                break  
    
    return sp.Add(*filtered_terms)

c1 = sp.Rational(1, 2)
c2 = sp.Rational(1, 6)
c3 = sp.Rational(1, 24)

# Hausdorff Expansion of Similarity Transformed Hamiltonian, assumming single and double excitations
Hbar = H + comm(H, A + B) + c1 * comm(comm(H, A + B), A + B) + c2 * comm(comm(comm(H, A + B), A + B), A + B) + c3 * comm(comm(comm(comm(H, A + B), A + B), A + B), A + B)
#Hbar2 = H + comm(H, A+ B + C + D) + c1 * comm(comm(H, A + B + C + D), A + B + C + D) + c2 * comm(comm(comm(H, A + B + C + D), A + B + C + D), A + B + C + D) + c3 * comm(comm(comm(comm(H, A + B + C + D ), A + B + C + D), A + B + C + D), A + B + C + D)
Hbar3 = H + comm(H,A) + c1 * comm(comm(H,A),A) + c2 * comm(comm(comm(H,A),A),A) + c3 * comm(comm(comm(comm(H,A),A),A),A)
print("Original expression:")
sp.pprint(Hbar, use_unicode=True)

# Simplify expressions
expand = simpl(Hbar)
expand2 = simpl(Hbar3)

print("Simplified")
sp.pprint(expand, use_unicode = True)
filtered = filter(expand, H, A, B)
filtered2 = filter(expand2, H, A, B)

# Print simplified and filtered expression
print("\n Conected unsimplified:")
sp.pprint(filtered, use_unicode=True)

A_com, B_com = sp.symbols('A B', commutative = True)
A_c, B_c, C_c, D_c, E_c, F_c = sp.symbols('A B C D E F', commutative = True)
simp = filtered.subs({A: A_com, B: B_com})
simp2 = filtered2.subs({A: A_c})
simp_exp = sp.simplify(simp)
simp_exp2 = sp.simplify(simp2)

print("\n Connected Hamiltonian (N)")
sp.pprint(simp_exp, use_unicode = True)

print("\n Connected CCSD-1-S1 Guess")
sp.pprint(simp_exp2, use_unicode = True)