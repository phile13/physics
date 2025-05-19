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
        const h = 6.62607015e-34;   // Planck constant (J·s)
        const c = 299792458;        // Speed of light (m/s)
        return h / (mass * c);
    }

    /**
     * Calculates the relativistic energy of a particle.
     * E = γ m c^2 where γ = 1 / sqrt(1 - v^2/c^2)
     * @param {number} mass - Rest mass in kg
     * @param {number} velocity - Particle velocity in m/s
     * @returns {number} Relativistic energy in Joules (J)
     */
    static RelativisticEnergy(mass, velocity) {
        const c = 299792458;
        const gamma = 1 / Math.sqrt(1 - (velocity * velocity) / (c * c));
        return gamma * mass * c * c;
    }

    /**
     * Calculates the de Broglie wavelength.
     * λ = h / p
     * @param {number} momentum - Particle momentum in kg·m/s
     * @returns {number} Wavelength in meters (m)
     */
    static DeBroglieWavelength(momentum) {
        const h = 6.62607015e-34;
        return h / momentum;
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
        const k = 8.9875517923e9; // Coulomb's constant (N·m²/C²)
        return k * q / (r * r);
    }

    /**
     * Calculates the magnetic field at center of a circular loop.
     * B = (μ₀ * I) / (2 * R)
     * @param {number} I - Current in Amperes (A)
     * @param {number} R - Radius of the loop in meters (m)
     * @returns {number} Magnetic field in Tesla (T)
     */
    static MagneticFieldLoop(I, R) {
        const mu0 = 4 * Math.PI * 1e-7; // Vacuum permeability (T·m/A)
        return (mu0 * I) / (2 * R);
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
     * @param {number} g - Acceleration due to gravity in m/s^2 (default 9.81)
     * @returns {number} Potential energy in Joules (J)
     */
    static GravitationalPotentialEnergy(mass, height, g = 9.81) {
        return mass * g * height;
    }

    /**
     * Calculates the period of a simple pendulum (small angles).
     * T = 2π * sqrt(L / g)
     * @param {number} length - Length of the pendulum in meters (m)
     * @param {number} g - Acceleration due to gravity in m/s^2 (default 9.81)
     * @returns {number} Period in seconds (s)
     */
    static PendulumPeriod(length, g = 9.81) {
        return 2 * Math.PI * Math.sqrt(length / g);
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
     * @param {number} g - Gravitational acceleration in m/s^2 (default 9.81)
     * @returns {number} Pressure at point 2 in Pascals (Pa)
     */
    static BernoulliPressure(P1, rho, v1, h1, v2, h2, g = 9.81) {
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
        const G = 6.67430e-11; // Gravitational constant (m^3·kg^-1·s^-2)
        return G * m1 * m2 / (r * r);
    }

    /**
     * Calculates the escape velocity from a spherical body.
     * v = sqrt(2GM / R)
     * @param {number} M - Mass of the body in kilograms (kg)
     * @param {number} R - Radius of the body in meters (m)
     * @returns {number} Escape velocity in meters per second (m/s)
     */
    static EscapeVelocity(M, R) {
        const G = 6.67430e-11;
        return Math.sqrt(2 * G * M / R);
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
        const G = 6.67430e-11;
        return -G * m1 * m2 / r;
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
        return 2 * Math.PI * frequency;
    }

    /**
     * Calculates the wave number.
     * k = 2π / λ
     * @param {number} wavelength - Wavelength in meters (m)
     * @returns {number} Wave number in radians per meter (rad/m)
     */
    static WaveNumber(wavelength) {
        return 2 * Math.PI / wavelength;
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
        const h = 6.62607015e-34;
        const c = 299792458;
        return (h * c) / wavelength;
    }

    /**
     * Calculates photon momentum.
     * p = h / λ
     * @param {number} wavelength - Wavelength in meters (m)
     * @returns {number} Momentum in kg·m/s
     */
    static PhotonMomentum(wavelength) {
        const h = 6.62607015e-34;
        return h / wavelength;
    }

    /**
     * Calculates the refractive index from speed of light in the medium.
     * n = c / v
     * @param {number} speedInMedium - Speed of light in the medium (m/s)
     * @returns {number} Refractive index (dimensionless)
     */
    static RefractiveIndex(speedInMedium) {
        const c = 299792458;
        return c / speedInMedium;
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
        const c = 299792458;
        return deltaT / Math.sqrt(1 - (velocity * velocity) / (c * c));
    }

    /**
     * Calculates length contraction.
     * L' = L * sqrt(1 - v^2/c^2)
     * @param {number} length - Proper length in meters (m)
     * @param {number} velocity - Relative velocity in m/s
     * @returns {number} Contracted length in meters (m)
     */
    static LengthContraction(length, velocity) {
        const c = 299792458;
        return length * Math.sqrt(1 - (velocity * velocity) / (c * c));
    }

    /**
     * Calculates relativistic momentum.
     * p = γ m v
     * @param {number} mass - Rest mass in kg
     * @param {number} velocity - Velocity in m/s
     * @returns {number} Relativistic momentum in kg·m/s
     */
    static RelativisticMomentum(mass, velocity) {
        const c = 299792458;
        const gamma = 1 / Math.sqrt(1 - (velocity * velocity) / (c * c));
        return gamma * mass * velocity;
    }
    
    /**
     * Calculates relativistic kinetic energy.
     * KE = (γ - 1) m c^2
     * @param {number} mass - Rest mass in kg
     * @param {number} velocity - Velocity in m/s
     * @returns {number} Kinetic energy in Joules (J)
     */
    static RelativisticKineticEnergy(mass, velocity) {
        const c = 299792458;
        const gamma = 1 / Math.sqrt(1 - (velocity * velocity) / (c * c));
        return (gamma - 1) * mass * c * c;
    }
}
