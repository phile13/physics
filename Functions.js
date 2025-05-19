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

class ClassicCollisionFormulas {
    /**
     * Calculates the final velocities after a 2D partially inelastic collision between two particles.
     * The coefficient of restitution e (0 ≤ e ≤ 1) determines elasticity.
     * @param {number} m1 - Mass of particle 1 (kg)
     * @param {number} m2 - Mass of particle 2 (kg)
     * @param {number} v1x - Initial x velocity of particle 1 (m/s)
     * @param {number} v1y - Initial y velocity of particle 1 (m/s)
     * @param {number} v2x - Initial x velocity of particle 2 (m/s)
     * @param {number} v2y - Initial y velocity of particle 2 (m/s)
     * @param {number} x1 - x position of particle 1 (m)
     * @param {number} y1 - y position of particle 1 (m)
     * @param {number} x2 - x position of particle 2 (m)
     * @param {number} y2 - y position of particle 2 (m)
     * @param {number} e  - Coefficient of restitution (0 for perfectly inelastic(stick together), 1 for perfectly elastic(bounce off each other))
     * @returns {object} {v1fx, v1fy, v2fx, v2fy} - Final velocities
     */
    static Collision2d(m1, m2, v1x, v1y, v2x, v2y, x1, y1, x2, y2, e) {
        // Normal vector between centers
        const dx = x2 - x1;
        const dy = y2 - y1;
        const distSq = dx*dx + dy*dy;
        if (distSq === 0) {
            // Overlapping particles, return unchanged
            return {v1fx: v1x, v1fy: v1y, v2fx: v2x, v2fy: v2y};
        }
        const dist = Math.sqrt(distSq);
        const nx = dx / dist;
        const ny = dy / dist;

        // Relative velocity
        const dvx = v1x - v2x;
        const dvy = v1y - v2y;

        // Velocity along the normal
        const relVel = dvx * nx + dvy * ny;

        // If particles are moving apart, don't collide
        if (relVel > 0) return {v1fx: v1x, v1fy: v1y, v2fx: v2x, v2fy: v2y};

        // Coefficient of restitution formula (see e.g. https://en.wikipedia.org/wiki/Elastic_collision)
        // v1n' = ( (m1 - e*m2)*v1n + (1+e)*m2*v2n ) / (m1 + m2)
        // v2n' = ( (m2 - e*m1)*v2n + (1+e)*m1*v1n ) / (m1 + m2)
        // We'll decompose velocities into normal and tangential components

        // Normal components
        const v1n = v1x * nx + v1y * ny;
        const v2n = v2x * nx + v2y * ny;
        // Tangential components (unchanged)
        const v1t = -v1x * ny + v1y * nx;
        const v2t = -v2x * ny + v2y * nx;

        // New normal velocities after collision
        const v1n_prime = ((m1 - e * m2) * v1n + (1 + e) * m2 * v2n) / (m1 + m2);
        const v2n_prime = ((m2 - e * m1) * v2n + (1 + e) * m1 * v1n) / (m1 + m2);

        // Convert back to x/y
        const v1fx = v1n_prime * nx - v1t * ny;
        const v1fy = v1n_prime * ny + v1t * nx;
        const v2fx = v2n_prime * nx - v2t * ny;
        const v2fy = v2n_prime * ny + v2t * nx;

        return {v1fx, v1fy, v2fx, v2fy};
    }

    /**
     * Calculates impulse delivered during a collision.
     * J = m * (v_after - v_before)
     * @param {number} mass - Mass in kg
     * @param {number} v_before - Velocity before collision (m/s)
     * @param {number} v_after - Velocity after collision (m/s)
     * @returns {number} Impulse in kg·m/s
     */
    static Impulse(mass, v_before, v_after) {
        return mass * (v_after - v_before);
    }

    /**
     * Calculates the kinetic energy lost in a collision.
     * ΔKE = KE_initial - KE_final
     * @param {number} m1 - Mass of particle 1 (kg)
     * @param {number} v1i - Initial velocity of particle 1 (m/s)
     * @param {number} m2 - Mass of particle 2 (kg)
     * @param {number} v2i - Initial velocity of particle 2 (m/s)
     * @param {number} v1f - Final velocity of particle 1 (m/s)
     * @param {number} v2f - Final velocity of particle 2 (m/s)
     * @returns {number} Kinetic energy lost (J)
     */
    static KineticEnergyLoss(m1, v1i, m2, v2i, v1f, v2f) {
        const KEi = 0.5 * m1 * v1i * v1i + 0.5 * m2 * v2i * v2i;
        const KEf = 0.5 * m1 * v1f * v1f + 0.5 * m2 * v2f * v2f;
        return KEi - KEf;
    }

    /**
     * Calculates the coefficient of restitution from measured velocities.
     * e = (v2f - v1f) / (v1i - v2i)
     * @param {number} v1i - Initial velocity of particle 1 (m/s)
     * @param {number} v2i - Initial velocity of particle 2 (m/s)
     * @param {number} v1f - Final velocity of particle 1 (m/s)
     * @param {number} v2f - Final velocity of particle 2 (m/s)
     * @returns {number} Coefficient of restitution (dimensionless)
     */
    static CoefficientOfRestitution(v1i, v2i, v1f, v2f) {
        return (v2f - v1f) / (v1i - v2i);
    }

    /**
     * Checks 2D momentum conservation before and after collision.
     * @param {number} m1 - Mass of particle 1 (kg)
     * @param {number} v1ix - Initial x velocity of particle 1 (m/s)
     * @param {number} v1iy - Initial y velocity of particle 1 (m/s)
     * @param {number} m2 - Mass of particle 2 (kg)
     * @param {number} v2ix - Initial x velocity of particle 2 (m/s)
     * @param {number} v2iy - Initial y velocity of particle 2 (m/s)
     * @param {number} m1f - Mass of particle 1 after (kg)
     * @param {number} v1fx - Final x velocity of particle 1 (m/s)
     * @param {number} v1fy - Final y velocity of particle 1 (m/s)
     * @param {number} m2f - Mass of particle 2 after (kg)
     * @param {number} v2fx - Final x velocity of particle 2 (m/s)
     * @param {number} v2fy - Final y velocity of particle 2 (m/s)
     * @returns {object} { px_i, py_i, px_f, py_f }
     */
    static MomentumConservation2D(
        m1, v1ix, v1iy, m2, v2ix, v2iy,
        m1f, v1fx, v1fy, m2f, v2fx, v2fy
    ) {
        const px_i = m1 * v1ix + m2 * v2ix;
        const py_i = m1 * v1iy + m2 * v2iy;
        const px_f = m1f * v1fx + m2f * v2fx;
        const py_f = m1f * v1fy + m2f * v2fy;
        return { px_i, py_i, px_f, py_f };
    }

    /**
     * Calculates the Mandelstam s variable for relativistic collisions.
     * s = (E1 + E2)^2 - (px1 + px2)^2 - (py1 + py2)^2 - (pz1 + pz2)^2
     * @param {number} E1 - Energy of particle 1 (J)
     * @param {number} px1 - x momentum of particle 1 (kg·m/s)
     * @param {number} py1 - y momentum of particle 1 (kg·m/s)
     * @param {number} pz1 - z momentum of particle 1 (kg·m/s)
     * @param {number} E2 - Energy of particle 2 (J)
     * @param {number} px2 - x momentum of particle 2 (kg·m/s)
     * @param {number} py2 - y momentum of particle 2 (kg·m/s)
     * @param {number} pz2 - z momentum of particle 2 (kg·m/s)
     * @returns {number} Mandelstam s (J^2)
     */
    static MandelstamS(E1, px1, py1, pz1, E2, px2, py2, pz2) {
        const E = E1 + E2;
        const px = px1 + px2;
        const py = py1 + py2;
        const pz = pz1 + pz2;
        return E*E - (px*px + py*py + pz*pz * Constants.SPEED_OF_LIGHT_SQUARED);
    }

    /**
     * Calculates the lab frame scattering angle from CM frame angle.
     * θ_lab = arctan( sinθ_cm / (γ + cosθ_cm) )
     * @param {number} theta_cm - Scattering angle in CM frame (radians)
     * @param {number} gamma - Dimensionless parameter (depends on masses/velocities)
     * @returns {number} θ_lab in radians
     */
    static LabScatteringAngle(theta_cm, gamma) {
        return Math.atan( Math.sin(theta_cm) / (gamma + Math.cos(theta_cm)) );
    }

    /**
     * Calculates the impact parameter for a collision.
     * b = |r0x * vy - r0y * vx| / sqrt(vx^2 + vy^2)
     * @param {number} r0x - Initial x position (m)
     * @param {number} r0y - Initial y position (m)
     * @param {number} vx - Velocity x (m/s)
     * @param {number} vy - Velocity y (m/s)
     * @returns {number} Impact parameter (m)
     */
    static ImpactParameter(r0x, r0y, vx, vy) {
        return Math.abs(r0x * vy - r0y * vx) / Math.sqrt(vx*vx + vy*vy);
    }

    /**
     * Calculates mean free path for collisions in a medium.
     * λ = 1 / (n * σ)
     * @param {number} numberDensity - Number density (particles/m^3)
     * @param {number} crossSection - Cross-section (m^2)
     * @returns {number} Mean free path (m)
     */
    static MeanFreePath(numberDensity, crossSection) {
        return 1 / (numberDensity * crossSection);
    }

    /**
     * Calculates the average time between collisions.
     * τ = λ / v
     * @param {number} meanFreePath - Mean free path (m)
     * @param {number} velocity - Particle velocity (m/s)
     * @returns {number} Time between collisions (s)
     */
    static CollisionTime(meanFreePath, velocity) {
        return meanFreePath / velocity;
    }
}

class QuantumScattering2D {
    /**
     * Simulate a 2D quantum elastic collision between two particles,
     * each with {mass, vx, vy}.
     * @param {object} p1 - {mass, vx, vy}
     * @param {object} p2 - {mass, vx, vy}
     * @returns {object} {particle1: {vx, vy}, particle2: {vx, vy}}
     */
    static Scatter(p1, p2) {
        const m1 = p1.mass, m2 = p2.mass;
        // 1. Compute momenta
        const mom1 = [m1 * p1.vx, m1 * p1.vy];
        const mom2 = [m2 * p2.vx, m2 * p2.vy];

        // 2. Total momentum and center-of-mass velocity
        const Ptot = [mom1[0] + mom2[0], mom1[1] + mom2[1]];
        const v_cm = [Ptot[0] / (m1 + m2), Ptot[1] / (m1 + m2)];

        // 3. Momenta in CM frame
        function subtract(a, b) { return [a[0] - b[0], a[1] - b[1]]; }
        const mom1_cm = subtract(mom1, [m1 * v_cm[0], m1 * v_cm[1]]);
        const mom2_cm = subtract(mom2, [m2 * v_cm[0], m2 * v_cm[1]]);

        // 4. Magnitude of outgoing momentum (conserved)
        const p_mag = Math.sqrt(mom1_cm[0] ** 2 + mom1_cm[1] ** 2);

        // 5. Random outgoing direction (isotropic in 2D)
        const theta = 2 * Math.PI * Math.random();
        const dir = [Math.cos(theta), Math.sin(theta)];

        // 6. New momenta in CM frame
        const mom1_cm_new = [p_mag * dir[0], p_mag * dir[1]];
        const mom2_cm_new = [-p_mag * dir[0], -p_mag * dir[1]];

        // 7. Convert momenta to velocities in CM frame
        const v1_cm_new = [mom1_cm_new[0] / m1, mom1_cm_new[1] / m1];
        const v2_cm_new = [mom2_cm_new[0] / m2, mom2_cm_new[1] / m2];

        // 8. Boost velocities back to lab frame
        const v1_lab = [v1_cm_new[0] + v_cm[0], v1_cm_new[1] + v_cm[1]];
        const v2_lab = [v2_cm_new[0] + v_cm[0], v2_cm_new[1] + v_cm[1]];

        return {
            particle1: {vx: v1_lab[0], vy: v1_lab[1]},
            particle2: {vx: v2_lab[0], vy: v2_lab[1]}
        };
    }
}

class QuantumRelativisticEvent2D {
    // --- VECTOR HELPERS ---
    static vecAdd(a, b) { return [a[0] + b[0], a[1] + b[1]]; }
    static vecSub(a, b) { return [a[0] - b[0], a[1] - b[1]]; }
    static vecScale(a, s) { return [a[0] * s, a[1] * s]; }
    static vecDot(a, b) { return a[0] * b[0] + a[1] * b[1]; }
    static vecNorm(a) { return Math.sqrt(a[0] * a[0] + a[1] * a[1]); }
    static vecUnit(a) { const n = this.vecNorm(a); return n === 0 ? [1, 0] : [a[0] / n, a[1] / n]; }

    // --- RELATIVISTIC HELPERS ---
    static momentum(m, v) {
        const vmag = this.vecNorm(v);
        return RelativisticFormulas.RelativisticMomentum(m, vmag);
    }
    static energy(m, v) {
        const vmag = this.vecNorm(v);
        return ParticleFormulas.RelativisticEnergy(m, vmag);
    }

    // --- LORENTZ BOOSTS ---
    static lorentzBoost(px, py, E, v_cm) {
        const v_cm_mag = this.vecNorm(v_cm);
        if (v_cm_mag === 0) return { px, py, E };
        const c = Constants.SPEED_OF_LIGHT;
        const beta = v_cm_mag / c;
        const gamma = 1 / Math.sqrt(1 - beta * beta);
        const n = this.vecUnit(v_cm);
        const p_par = px * n[0] + py * n[1];
        const p_perp = [-n[1] * px + n[0] * py];
        const p_par_prime = gamma * (p_par - beta * E / c);
        const E_prime = gamma * (E - beta * c * p_par);
        const px_prime = p_par_prime * n[0] - p_perp[0] * n[1];
        const py_prime = p_par_prime * n[1] + p_perp[0] * n[0];
        return { px: px_prime, py: py_prime, E: E_prime };
    }
    static inverseLorentzBoost(px, py, E, v_cm) {
        const v_cm_mag = this.vecNorm(v_cm);
        if (v_cm_mag === 0) return { px, py, E };
        const c = Constants.SPEED_OF_LIGHT;
        const beta = v_cm_mag / c;
        const gamma = 1 / Math.sqrt(1 - beta * beta);
        const n = this.vecUnit(v_cm);
        const p_par = px * n[0] + py * n[1];
        const p_perp = [-n[1] * px + n[0] * py];
        const p_par_prime = gamma * (p_par + beta * E / c);
        const E_prime = gamma * (E + beta * c * p_par);
        const px_prime = p_par_prime * n[0] - p_perp[0] * n[1];
        const py_prime = p_par_prime * n[1] + p_perp[0] * n[0];
        return { px: px_prime, py: py_prime, E: E_prime };
    }
    static velocityFromMomentum(px, py, mass) {
        const p = Math.sqrt(px ** 2 + py ** 2);
        const c = Constants.SPEED_OF_LIGHT;
        const gamma = Math.sqrt(1 + (p / (mass * c)) ** 2);
        const v = p / (mass * gamma);
        const norm = Math.sqrt(px ** 2 + py ** 2);
        return norm === 0 ? [0, 0] : [v * px / norm, v * py / norm];
    }

    // --- PHASE SPACE GENERATOR (RAMBO-LIKE) ---
    static generatePhaseSpace2D(N, ECM, masses) {
        // 1. Generate random directions and energies for massless particles in CM
        let q = [];
        let sum_px = 0, sum_py = 0, sum_E = 0;
        for (let i = 0; i < N; ++i) {
            const theta = 2 * Math.PI * Math.random();
            const xi = -Math.log(Math.random());
            const px = xi * Math.cos(theta);
            const py = xi * Math.sin(theta);
            const E = xi;
            q.push({ px, py, E });
            sum_px += px;
            sum_py += py;
            sum_E += E;
        }
        // 2. Center momenta
        for (let i = 0; i < N; ++i) {
            q[i].px -= sum_px / N;
            q[i].py -= sum_py / N;
        }
        // 3. Rescale energies
        sum_E = 0;
        for (let i = 0; i < N; ++i) {
            q[i].E = Math.sqrt(q[i].px ** 2 + q[i].py ** 2);
            sum_E += q[i].E;
        }
        const scale = ECM / sum_E;
        for (let i = 0; i < N; ++i) {
            q[i].px *= scale;
            q[i].py *= scale;
            q[i].E *= scale;
        }
        // 4. Assign masses and correct momenta/energies
        const c2 = Constants.SPEED_OF_LIGHT_SQUARED;
        for (let i = 0; i < N; ++i) {
            const m = masses[i];
            const p_mag = Math.sqrt(q[i].px ** 2 + q[i].py ** 2);
            const E = Math.sqrt(p_mag ** 2 * c2 + m ** 2 * c2 * c2);
            const p_desired = Math.sqrt((E ** 2 - m ** 2 * c2 * c2) / c2);
            q[i].px = q[i].px * (p_desired / p_mag);
            q[i].py = q[i].py * (p_desired / p_mag);
            q[i].E = E;
        }
        return q;
    }

    // --- MAIN EVENT GENERATOR ---
    static selectAndScatter(particleA, particleB, outcomes) {
        // 1. Weighted random selection
        const totalProb = outcomes.reduce((sum, o) => sum + o.prob, 0);
        let r = Math.random() * totalProb;
        let selected = outcomes[0];
        for (const o of outcomes) {
            if (r < o.prob) { selected = o; break; }
            r -= o.prob;
        }

        // 2. Compute total momentum and energy in LAB
        const mA = particleA.mass, mB = particleB.mass;
        const vA = [particleA.vx, particleA.vy], vB = [particleB.vx, particleB.vy];
        const pA_mag = this.momentum(mA, vA);
        const pB_mag = this.momentum(mB, vB);
        const pA_dir = this.vecUnit(vA);
        const pB_dir = this.vecUnit(vB);
        const pA_vec = this.vecScale(pA_dir, pA_mag);
        const pB_vec = this.vecScale(pB_dir, pB_mag);
        const pTot = this.vecAdd(pA_vec, pB_vec);
        const eA = this.energy(mA, vA);
        const eB = this.energy(mB, vB);
        const ETot = eA + eB;

        // 3. Center-of-mass velocity (relativistic)
        const c2 = Constants.SPEED_OF_LIGHT_SQUARED;
        const v_cm = [pTot[0] * c2 / ETot, pTot[1] * c2 / ETot];

        // 4. Boost incoming particles to CM frame
        const a_cm = this.lorentzBoost(pA_vec[0], pA_vec[1], eA, v_cm);
        const b_cm = this.lorentzBoost(pB_vec[0], pB_vec[1], eB, v_cm);
        const ECM = a_cm.E + b_cm.E;

        // 5. Generate final state
        let newParticles = [];
        const N = selected.particles.length;
        if (N === 2) {
            // --- 2-body final state (as before) ---
            const m1 = selected.particles[0].mass, m2 = selected.particles[1].mass;
            const s = ECM * ECM / c2 / c2;
            const term1 = s - Math.pow((m1 + m2), 2);
            const term2 = s - Math.pow((m1 - m2), 2);
            const p_star = 0.5 * Math.sqrt(term1 * term2) * Constants.SPEED_OF_LIGHT;
            const theta = 2 * Math.PI * Math.random();
            const dir = [Math.cos(theta), Math.sin(theta)];
            const p1_cm = [p_star * dir[0], p_star * dir[1]];
            const p2_cm = [-p_star * dir[0], -p_star * dir[1]];
            const E1_cm = Math.sqrt((p_star * p_star * c2) + (m1 * c2) * (m1 * c2));
            const E2_cm = Math.sqrt((p_star * p_star * c2) + (m2 * c2) * (m2 * c2));
            const p1_lab = this.inverseLorentzBoost(p1_cm[0], p1_cm[1], E1_cm, v_cm);
            const p2_lab = this.inverseLorentzBoost(p2_cm[0], p2_cm[1], E2_cm, v_cm);
            const v1 = this.velocityFromMomentum(p1_lab.px, p1_lab.py, m1);
            const v2 = this.velocityFromMomentum(p2_lab.px, p2_lab.py, m2);
            newParticles = [
                { mass: m1, vx: v1[0], vy: v1[1] },
                { mass: m2, vx: v2[0], vy: v2[1] }
            ];
        } else {
            // --- N-body phase space (RAMBO) ---
            const masses = selected.particles.map(p => p.mass);
            const q = this.generatePhaseSpace2D(N, ECM, masses);
            for (let i = 0; i < N; ++i) {
                const m = masses[i];
                const lab = this.inverseLorentzBoost(q[i].px, q[i].py, q[i].E, v_cm);
                const v = this.velocityFromMomentum(lab.px, lab.py, m);
                newParticles.push({ mass: m, vx: v[0], vy: v[1] });
            }
        }
        return { selectedOutcome: selected, newParticles };
    }
}



