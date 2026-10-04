def area_mach_residual(M, target_AR, gamma):
    # The standard Area-Mach relation
    term1 = 2 / (gamma + 1)
    term2 = 1 + (gamma - 1) / 2 * M**2
    exponent = (gamma + 1) / (2 * (gamma - 1))

    calculated_AR = (1/M) * (term1 * term2)**exponent
    return calculated_AR - target_AR
