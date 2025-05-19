class Constants {
    // 2 PI
    static TWO_PI = 2 * Math.PI
    
    // Planck constant (Joule seconds)
    static PLANCK = 6.62607015e-34;

    // Reduced Planck constant (ħ = h / 2π) (Joule seconds)
    static HBAR = 6.62607015e-34 / (Constants.TWO_PI);

    // Speed of light in vacuum (meters per second)
    static SPEED_OF_LIGHT = 299792458;

    // Speed of light squared (m^2/s^2)
    static SPEED_OF_LIGHT_SQUARED = 299792458 ** 2;

    // Speed of light cubed (m^3/s^3)
    static SPEED_OF_LIGHT_CUBED = 299792458 ** 3;

    // Gravitational constant (m^3 kg^-1 s^-2)
    static GRAVITATIONAL = 6.67430e-11;

    // Coulomb's constant (N·m²/C²)
    static COULOMB = 8.9875517923e9;

    // Vacuum permeability (Tesla meter per Ampere)
    static VACUUM_PERMEABILITY = 4 * Math.PI * 1e-7;

    // Standard acceleration due to gravity (m/s^2)
    static GRAVITY = 9.80665;

    // Elementary charge (Coulombs)
    static ELECTRON_CHARGE = 1.602176634e-19;

    // Electron rest mass (kg)
    static ELECTRON_MASS = 9.10938356e-31;

    // Boltzmann constant (Joule per Kelvin)
    static BOLTZMANN = 1.380649e-23;

    // Ideal gas constant (Joule per mole per Kelvin)
    static IDEAL_GAS = 8.314462618;

    // Avogadro's number (particles per mole)
    static AVOGADRO = 6.02214076e23;

    // Stefan-Boltzmann constant (W m^-2 K^-4)
    static STEFAN_BOLTZMANN = 5.670374419e-8;

    // Proton rest mass (kg)
    static PROTON_MASS = 1.67262192369e-27;
    
    // Neutron rest mass (kg)
    static NEUTRON_MASS = 1.67492749804e-27;
    
    // Atomic mass unit (kg)
    static AMU = 1.66053906660e-27;
}

class NuclearFormulas {
    /**
     * Calculates missing mass in a nuclear reaction (from search result [2])
     * m² = (E_initial - E_final_known)² - |p_initial - p_final_known|²
     * @param {number} E_initial - Total initial energy (J)
     * @param {number} E_final_known - Sum of known final energies (J)
     * @param {number} p_initial - Initial momentum magnitude (kg·m/s)
     * @param {number} p_final_known - Sum of known final momenta magnitudes (kg·m/s)
     * @returns {number} Missing mass (kg)
     */
    static MissingMass(E_initial, E_final_known, p_initial, p_final_known) {
        const energyDiff = E_initial - E_final_known;
        const momentumDiff = p_initial - p_final_known;
        return Math.sqrt(energyDiff**2 - (momentumDiff * Constants.SPEED_OF_LIGHT)**2) 
               / Constants.SPEED_OF_LIGHT_SQUARED;
    }

    /**
     * Calculates decay constant from half-life (from search result [6])
     * λ = ln(2) / T½
     * @param {number} halfLife - Half-life in seconds
     * @returns {number} Decay constant (s⁻¹)
     */
    static DecayConstant(halfLife) {
        return Math.LN2 / halfLife;
    }

    /**
     * Calculates remaining nuclei after time t (from search result [6])
     * N = N₀e^(-λt)
     * @param {number} N0 - Initial quantity
     * @param {number} lambda - Decay constant (s⁻¹)
     * @param {number} t - Time elapsed (s)
     * @returns {number} Remaining nuclei
     */
    static RemainingNuclei(N0, lambda, t) {
        return N0 * Math.exp(-lambda * t);
    }

    /**
     * Alpha decay equation (from search result [5])
     * Returns daughter nucleus and alpha particle masses
     * @param {number} Z - Atomic number of parent
     * @param {number} A - Mass number of parent
     * @returns {Object} {daughterZ, daughterA, alphaMass}
     */
    static AlphaDecay(Z, A) {
        return {
            daughterZ: Z - 2,
            daughterA: A - 4,
            alphaMass: 4 * Constants.PROTON_MASS + 2 * Constants.ELECTRON_MASS
        };
    }
    
    /**
     * Binding energy formula (Bethe-Weizsäcker approximation)
     * B = aV*A - aS*A^(2/3) - aC*Z²/A^(1/3) - aA*(A-2Z)²/A + δ(A,Z)
     * @param {number} A - Mass number
     * @param {number} Z - Atomic number
     * @returns {number} Binding energy in MeV
     */
    static BindingEnergy(A, Z) {
        const aV = 15.5, aS = 16.8, aC = 0.715, aA = 23.0; // Coefficients in MeV
        const surface = aS * Math.pow(A, 2/3);
        const coulomb = aC * Z**2 / Math.pow(A, 1/3);
        const asymmetry = aA * (A - 2*Z)**2 / A;
        const pairing = (A%2 === 0 && Z%2 === 0) ? 12/Math.sqrt(A) : 0;
        return aV*A - surface - coulomb - asymmetry + pairing;
    }

    /**
     * Neutron moderation equation (from search result [3])
     * ξ = ln(E_initial/E_final) - Average logarithmic energy loss per collision
     * @param {number} A - Mass number of moderator nucleus
     * @returns {number} Average logarithmic energy decrement
     */
    static NeutronModeration(A) {
        return 1 + ((A-1)**2/(2*A)) * Math.log((A-1)/(A+1));
    }

    /**
     * Q-value calculation (from search result [3])
     * Q = (Σm_initial - Σm_final)c²
     * @param {number[]} initialMasses - Array of initial particle masses (kg)
     * @param {number[]} finalMasses - Array of final particle masses (kg)
     * @returns {number} Q-value in joules
     */
    static QValue(initialMasses, finalMasses) {
        const sumInitial = initialMasses.reduce((a,b) => a + b, 0);
        const sumFinal = finalMasses.reduce((a,b) => a + b, 0);
        return (sumInitial - sumFinal) * Constants.SPEED_OF_LIGHT_SQUARED;
    }
}

class ParticleFormulas {
    /**
     * Calculates electric charge from isospin and hypercharge using the Gell-Mann–Nishijima formula.
     * Q = I3 + Y/2
     * @param {number} I3 - Third component of isospin
     * @param {number} Y - Hypercharge
     * @returns {number} Electric charge in units of e
     */
    static ElectricCharge(I3, Y) {
        return I3 + Y / 2;
    }

    /**
     * Calculates the Compton wavelength of a particle.
     * λ = h / (m * c)
     * @param {number} mass - Particle rest mass in kilograms (kg)
     * @returns {number} Compton wavelength in meters (m)
     */
    static ComptonWavelength(mass) {
        return Constants.PLANCK / (mass * Constants.SPEED_OF_LIGHT);
    }

    /**
     * Calculates the relativistic energy of a particle.
     * E = γ m c^2 where γ = 1 / sqrt(1 - v^2/c^2)
     * @param {number} mass - Rest mass in kg
     * @param {number} velocity - Particle velocity in m/s
     * @returns {number} Relativistic energy in Joules (J)
     */
    static RelativisticEnergy(mass, velocity) {
        return mass * Constants.SPEED_OF_LIGHT_SQUARED / Math.sqrt(1 - (velocity * velocity) / Constants.SPEED_OF_LIGHT_SQUARED);
    }

    /**
     * Calculates the de Broglie wavelength.
     * λ = h / p
     * @param {number} momentum - Particle momentum in kg·m/s
     * @returns {number} Wavelength in meters (m)
     */
    static DeBroglieWavelength(momentum) {
        return Constants.PLANCK / momentum;
    }
    
    /**
     * Calculates beta decay energy (from search result [6])
     * Q = (m_parent - m_daughter)c²
     * @param {number} massParent - Parent nucleus mass (kg)
     * @param {number} massDaughter - Daughter nucleus mass (kg)
     * @returns {number} Maximum kinetic energy of products (J)
     */
    static BetaDecayEnergy(massParent, massDaughter) {
        return (massParent - massDaughter) * Constants.SPEED_OF_LIGHT_SQUARED;
    }

    /**
     * Calculates center-of-mass energy (from search result [4])
     * √s = √(E₁ + E₂)² - |p₁c + p₂c|²
     * @param {number} E1 - Energy of particle 1 (J)
     * @param {number} p1 - Momentum of particle 1 (kg·m/s)
     * @param {number} E2 - Energy of particle 2 (J)
     * @param {number} p2 - Momentum of particle 2 (kg·m/s)
     * @returns {number} Center-of-mass energy (J)
     */
    static CenterOfMassEnergy(E1, p1, E2, p2) {
        const energySum = E1 + E2;
        const momentumSum = p1 + p2;
        return Math.sqrt(energySum**2 - (momentumSum * Constants.SPEED_OF_LIGHT)**2);
    }

    /**
     * Mandelstam variables (from search result [2])
     * s + t + u = m1² + m2² + m3² + m4²
     * @param {number} s - Mandelstam s
     * @param {number} t - Mandelstam t
     * @param {number} u - Mandelstam u
     * @param {number[]} masses - Array of four particle masses [m1,m2,m3,m4]
     * @returns {boolean} Conservation check
     */
    static MandelstamCheck(s, t, u, masses) {
        const massSum = masses.reduce((sum,m) => sum + m**2, 0);
        return Math.abs(s + t + u - massSum) < 1e-6;
    }

    /**
     * Breit-Wigner resonance formula (from search result [8])
     * σ(E) = σ_max * (Γ²/4) / [(E-M)² + Γ²/4]
     * @param {number} E - Center-of-mass energy
     * @param {number} M - Resonance mass
     * @param {number} Γ - Resonance width
     * @param {number} σ_max - Peak cross-section
     * @returns {number} Cross-section
     */
    static BreitWigner(E, M, Γ, σ_max) {
        return σ_max * (Γ**2/4) / ((E - M)**2 + Γ**2/4);
    }
}

class ElectromagnetismFormulas {
    /**
     * Calculates the electric field due to a point charge.
     * E = k * q / r^2
     * @param {number} q - Charge in Coulombs (C)
     * @param {number} r - Distance from charge in meters (m)
     * @returns {number} Electric field magnitude in N/C
     */
    static ElectricFieldPointCharge(q, r) {
        return Constants.COULOMB * q / (r * r);
    }

    /**
     * Calculates the magnetic field at center of a circular loop.
     * B = (μ₀ * I) / (2 * R)
     * @param {number} I - Current in Amperes (A)
     * @param {number} R - Radius of the loop in meters (m)
     * @returns {number} Magnetic field in Tesla (T)
     */
    static MagneticFieldLoop(I, R) {
        return (Constants.VACUUM_PERMEABILITY * I) / (2 * R);
    }

    /**
     * Calculates the force on a charge moving in a magnetic field.
     * F = q * v * B * sin(θ)
     * @param {number} q - Charge in Coulombs (C)
     * @param {number} v - Velocity in m/s
     * @param {number} B - Magnetic field in Tesla (T)
     * @param {number} theta - Angle between velocity and magnetic field in radians
     * @returns {number} Force in Newtons (N)
     */
    static LorentzForceMagnetic(q, v, B, theta) {
        return q * v * B * Math.sin(theta);
    }

    /**
     * Calculates the energy stored in a capacitor.
     * U = 1/2 C V^2
     * @param {number} C - Capacitance in Farads (F)
     * @param {number} V - Voltage in Volts (V)
     * @returns {number} Energy in Joules (J)
     */
    static CapacitorEnergy(C, V) {
        return 0.5 * C * V * V;
    }

    /**
     * Maxwell stress tensor (differential form)
     * T_ij = ε₀(E_iE_j - 0.5δ_ijE²) + (1/μ₀)(B_iB_j - 0.5δ_ijB²)
     * @param {number[]} E - Electric field vector [Ex,Ey,Ez]
     * @param {number[]} B - Magnetic field vector [Bx,By,Bz]
     * @returns {number[][]} 3x3 stress tensor
     */
    static MaxwellStressTensor(E, B) {
        const ε0 = Constants.VACUUM_PERMITTIVITY;
        const μ0 = Constants.VACUUM_PERMEABILITY;
        const E_sq = E.reduce((sum, e) => sum + e**2, 0);
        const B_sq = B.reduce((sum, b) => sum + b**2, 0);
        
        return Array.from({length:3}, (_,i) =>
            Array.from({length:3}, (_,j) => 
                ε0*(E[i]*E[j] - 0.5*(i===j)*E_sq) + 
                (1/μ0)*(B[i]*B[j] - 0.5*(i===j)*B_sq)
            )
        );
    }
}


class ClassicalMechanicsFormulas {
    /**
     * Calculates the kinetic energy of an object.
     * KE = 1/2 m v^2
     * @param {number} mass - Mass in kilograms (kg)
     * @param {number} velocity - Velocity in meters per second (m/s)
     * @returns {number} Kinetic energy in Joules (J)
     */
    static KineticEnergy(mass, velocity) {
        return 0.5 * mass * velocity * velocity;
    }

    /**
     * Calculates the gravitational potential energy near Earth's surface.
     * U = m g h
     * @param {number} mass - Mass in kilograms (kg)
     * @param {number} height - Height above reference in meters (m)
     * @param {number} g - Acceleration due to gravity in m/s^2 (default 9.80665)
     * @returns {number} Potential energy in Joules (J)
     */
    static GravitationalPotentialEnergy(mass, height, g = 9.80665) {
        return mass * g * height;
    }

    /**
     * Calculates the period of a simple pendulum (small angles).
     * T = 2π * sqrt(L / g)
     * @param {number} length - Length of the pendulum in meters (m)
     * @param {number} g - Acceleration due to gravity in m/s^2 (default 9.80665)
     * @returns {number} Period in seconds (s)
     */
    static PendulumPeriod(length, g = 9.80665) {
        return Constants.TWO_PI * Math.sqrt(length / g);
    }

    /**
     * Calculates the force using Hooke's Law.
     * F = -k x
     * @param {number} k - Spring constant in N/m
     * @param {number} displacement - Displacement from equilibrium in meters (m)
     * @returns {number} Restoring force in Newtons (N)
     */
    static HookesLawForce(k, displacement) {
        return -k * displacement;
    }
}

class FluidMechanicsFormulas {
    /**
     * Calculates the pressure difference using Bernoulli's equation.
     * P1 + 1/2 ρ v1^2 + ρ g h1 = P2 + 1/2 ρ v2^2 + ρ g h2
     * Returns pressure at point 2 assuming P1, v1, h1, v2, h2 known.
     * @param {number} P1 - Pressure at point 1 in Pascals (Pa)
     * @param {number} rho - Fluid density in kg/m^3
     * @param {number} v1 - Velocity at point 1 in m/s
     * @param {number} h1 - Height at point 1 in meters (m)
     * @param {number} v2 - Velocity at point 2 in m/s
     * @param {number} h2 - Height at point 2 in meters (m)
     * @param {number} g - Gravitational acceleration in m/s^2 (default 9.80665)
     * @returns {number} Pressure at point 2 in Pascals (Pa)
     */
    static BernoulliPressure(P1, rho, v1, h1, v2, h2, g = 9.80665) {
        return P1 + 0.5 * rho * v1 * v1 + rho * g * h1 - 0.5 * rho * v2 * v2 - rho * g * h2;
    }

    /**
     * Calculates flow rate using the continuity equation.
     * Q = A1 * v1 = A2 * v2
     * Returns velocity at point 2.
     * @param {number} A1 - Cross-sectional area at point 1 in m^2
     * @param {number} v1 - Velocity at point 1 in m/s
     * @param {number} A2 - Cross-sectional area at point 2 in m^2
     * @returns {number} Velocity at point 2 in m/s
     */
    static ContinuityVelocity(A1, v1, A2) {
        return (A1 * v1) / A2;
    }

    /**
     * Calculates the drag force on an object moving through a fluid.
     * Fd = 1/2 * ρ * v^2 * Cd * A
     * @param {number} rho - Fluid density in kg/m^3
     * @param {number} v - Velocity relative to fluid in m/s
     * @param {number} Cd - Drag coefficient (dimensionless)
     * @param {number} A - Cross-sectional area in m^2
     * @returns {number} Drag force in Newtons (N)
     */
    static DragForce(rho, v, Cd, A) {
        return 0.5 * rho * v * v * Cd * A;
    }
}

class GravitationFormulas {
    /**
     * Calculates the gravitational force between two masses.
     * F = G * m1 * m2 / r^2
     * @param {number} m1 - Mass 1 in kilograms (kg)
     * @param {number} m2 - Mass 2 in kilograms (kg)
     * @param {number} r - Distance between centers in meters (m)
     * @returns {number} Gravitational force in Newtons (N)
     */
    static GravitationalForce(m1, m2, r) {
        return Constants.GRAVITATIONAL * m1 * m2 / (r * r);
    }

    /**
     * Calculates the escape velocity from a spherical body.
     * v = sqrt(2GM / R)
     * @param {number} M - Mass of the body in kilograms (kg)
     * @param {number} R - Radius of the body in meters (m)
     * @returns {number} Escape velocity in meters per second (m/s)
     */
    static EscapeVelocity(M, R) {
        return Math.sqrt(2 * Constants.GRAVITATIONAL * M / R);
    }

    /**
     * Calculates gravitational potential energy between two masses.
     * U = - G * m1 * m2 / r
     * @param {number} m1 - Mass 1 in kg
     * @param {number} m2 - Mass 2 in kg
     * @param {number} r - Distance between masses in meters (m)
     * @returns {number} Potential energy in Joules (J)
     */
    static GravitationalPotentialEnergy(m1, m2, r) {
        return -Constants.GRAVITATIONAL * m1 * m2 / r;
    }
}

class WaveTheoryFormulas {
    /**
     * Calculates the wave speed.
     * v = f * λ
     * @param {number} frequency - Frequency in Hertz (Hz)
     * @param {number} wavelength - Wavelength in meters (m)
     * @returns {number} Wave speed in meters per second (m/s)
     */
    static WaveSpeed(frequency, wavelength) {
        return frequency * wavelength;
    }

    /**
     * Calculates the angular frequency.
     * ω = 2π f
     * @param {number} frequency - Frequency in Hertz (Hz)
     * @returns {number} Angular frequency in radians per second (rad/s)
     */
    static AngularFrequency(frequency) {
        return Constants.TWO_PI * frequency;
    }

    /**
     * Calculates the wave number.
     * k = 2π / λ
     * @param {number} wavelength - Wavelength in meters (m)
     * @returns {number} Wave number in radians per meter (rad/m)
     */
    static WaveNumber(wavelength) {
        return Constants.TWO_PI / wavelength;
    }

    /**
     * Calculates the intensity of a wave proportional to the square of amplitude.
     * I ∝ A^2
     * @param {number} amplitude - Wave amplitude (unit depends on context)
     * @returns {number} Intensity (arbitrary units)
     */
    static Intensity(amplitude) {
        return amplitude * amplitude;
    }
}

class PhotonicsFormulas {
    /**
     * Calculates photon energy from wavelength.
     * E = hc / λ
     * @param {number} wavelength - Wavelength in meters (m)
     * @returns {number} Photon energy in Joules (J)
     */
    static PhotonEnergy(wavelength) {
        return (Constants.PLANCK * Constants.SPEED_OF_LIGHT) / wavelength;
    }

    /**
     * Calculates photon momentum.
     * p = h / λ
     * @param {number} wavelength - Wavelength in meters (m)
     * @returns {number} Momentum in kg·m/s
     */
    static PhotonMomentum(wavelength) {
        return Constants.PLANCK / wavelength;
    }

    /**
     * Calculates the refractive index from speed of light in the medium.
     * n = c / v
     * @param {number} speedInMedium - Speed of light in the medium (m/s)
     * @returns {number} Refractive index (dimensionless)
     */
    static RefractiveIndex(speedInMedium) {
        return Constants.SPEED_OF_LIGHT / speedInMedium;
    }
}

class ThermodynamicsFormulas {
    /**
     * Entropy change for ideal gas
     * ΔS = nRln(V2/V1) + nCvln(T2/T1)
     * @param {number} n - Moles of gas
     * @param {number} V1 - Initial volume
     * @param {number} V2 - Final volume
     * @param {number} T1 - Initial temperature
     * @param {number} T2 - Final temperature
     * @param {number} Cv - Molar heat capacity at const volume
     * @returns {number} Entropy change (J/K)
     */
    static EntropyChange(n, V1, V2, T1, T2, Cv) {
        return n*Constants.IDEAL_GAS*Math.log(V2/V1) + n*Cv*Math.log(T2/T1);
    }

    /**
     * Enthalpy of formation
     * ΔH°_f = ΣΔH°_products - ΣΔH°_reactants
     * @param {number[]} productEnthalpies - Array of product ΔH° values
     * @param {number[]} reactantEnthalpies - Array of reactant ΔH° values
     * @returns {number} Enthalpy change (kJ/mol)
     */
    static FormationEnthalpy(productEnthalpies, reactantEnthalpies) {
        const sumProducts = productEnthalpies.reduce((a,b) => a + b, 0);
        const sumReactants = reactantEnthalpies.reduce((a,b) => a + b, 0);
        return sumProducts - sumReactants;
    }
}


class RelativisticFormulas {
    /**
     * Calculates time dilation.
     * Δt' = Δt / sqrt(1 - v^2/c^2)
     * @param {number} deltaT - Proper time interval in seconds (s)
     * @param {number} velocity - Relative velocity in m/s
     * @returns {number} Dilated time interval in seconds (s)
     */
    static TimeDilation(deltaT, velocity) {
        return deltaT * RelativisticFormulas.Gamma(velocity);
    }

    /**
     * Calculates length contraction.
     * L' = L * sqrt(1 - v^2/c^2)
     * @param {number} length - Proper length in meters (m)
     * @param {number} velocity - Relative velocity in m/s
     * @returns {number} Contracted length in meters (m)
     */
    static LengthContraction(length, velocity) {
        return length * RelativisticFormulas.GammaReciprocal(velocity);
    }
    
    /**
     * Calculates relativistic momentum.
     * γ = (1 - (v^2 /c^2))^-1/2
     * @param {number} velocity - Velocity in m/s
     * @returns {number} Relativistic momentum in kg·m/s
     */
    static GammaReciprocal(velocity) {
        return = Math.sqrt(1 - (velocity * velocity) / Constants.SPEED_OF_LIGHT_SQUARED);
    }
    
    /**
     * Calculates relativistic gamma function.
     * γ = (1 - (v^2 /c^2))^-1/2
     * @param {number} velocity - Velocity in m/s
     * @returns {number} Relativistic momentum in kg·m/s
     */
    static Gamma(velocity) {
        return = 1 / Math.sqrt(1 - (velocity * velocity) / Constants.SPEED_OF_LIGHT_SQUARED);
    }
    
    /**
     * Calculates relativistic momentum.
     * p = γ m v
     * @param {number} mass - Rest mass in kg
     * @param {number} velocity - Velocity in m/s
     * @returns {number} Relativistic momentum in kg·m/s
     */
    static RelativisticMomentum(mass, velocity) {
        return RelativisticFormulas.Gamma(velocity) * mass * velocity;
    }
    
    /**
     * Calculates relativistic kinetic energy.
     * KE = (γ - 1) m c^2
     * @param {number} mass - Rest mass in kg
     * @param {number} velocity - Velocity in m/s
     * @returns {number} Kinetic energy in Joules (J)
     */
    static RelativisticKineticEnergy(mass, velocity) {
        return (RelativisticFormulas.Gamma(velocity) - 1) * mass * Constants.SPEED_OF_LIGHT_SQUARED;
    }

    /**
     * Calculates two-body decay momentum (from search result [5])
     * p = √[(E² - (m₁ + m₂)²)(E² - (m₁ - m₂)²)] / (2E)
     * @param {number} E_total - Total energy of parent particle (J)
     * @param {number} m1 - Mass of first decay product (kg)
     * @param {number} m2 - Mass of second decay product (kg)
     * @returns {number} Momentum magnitude (kg·m/s)
     */
    static TwoBodyDecayMomentum(E_total, m1, m2) {
        const term1 = E_total**2/Constants.SPEED_OF_LIGHT_SQUARED - (m1 + m2)**2;
        const term2 = E_total**2/Constants.SPEED_OF_LIGHT_SQUARED - (m1 - m2)**2;
        return Math.sqrt(term1 * term2) / (2 * E_total/Constants.SPEED_OF_LIGHT_SQUARED);
    }

    /**
     * Calculates total energy from momentum (from search result [3])
     * E² = p²c² + m²c⁴
     * @param {number} momentum - Particle momentum (kg·m/s)
     * @param {number} mass - Rest mass (kg)
     * @returns {number} Total energy (J)
     */
    static EnergyFromMomentum(momentum, mass) {
        return Math.sqrt(
            (momentum * Constants.SPEED_OF_LIGHT)**2 + 
            (mass * Constants.SPEED_OF_LIGHT_SQUARED)**2
        );
    }

    /**
     * Rapidity calculation (from search result [11])
     * y = 0.5 * ln[(E+p)/(E-p)]
     * @param {number} E - Total energy (GeV)
     * @param {number} p - Momentum (GeV/c)
     * @returns {number} Rapidity
     */
    static Rapidity(E, p) {
        return 0.5 * Math.log((E + p)/(E - p));
    }

    /**
     * Invariant mass calculation (from search result [2])
     * M² = (ΣE)² - |Σp|²
     * @param {number[]} energies - Array of particle energies (GeV)
     * @param {number[][]} momenta - Array of 3-momentum vectors [px,py,pz] (GeV/c)
     * @returns {number} Invariant mass (GeV/c²)
     */
    static InvariantMass(energies, momenta) {
        const sumE = energies.reduce((a,b) => a + b, 0);
        const sumP = [0,0,0];
        momenta.forEach(p => {
            sumP[0] += p[0];
            sumP[1] += p[1];
            sumP[2] += p[2];
        });
        const pSq = sumP.reduce((a,b) => a + b**2, 0);
        return Math.sqrt(sumE**2 - pSq);
    }
}


