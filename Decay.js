class Decay {
    static CheckForDecay(particle){
        return (particle.age > Constants.ParticleInfo[particle.name].lifespan);
    }
    
    static HandleDecay(particle){
        if (Decay.CheckForDecay(particle)) {
            let combos = Decay.DetermineCombinations(particle);
            let sorted_combos = Object.values(combos).sort((a, b) => a.pct - b.pct);
    
            if (sorted_combos.length > 0) {
                // Calculate total probability sum for weighted selection
                let totalPct = sorted_combos.reduce((sum, combo) => sum + combo.pct, 0);
                
                // Generate a random number between 0 and totalPct
                let randomPick = Math.random() * totalPct;
    
                // Find the decay mode based on the random selection
                let cumulativePct = 0;
                let selectedCombo = null;
    
                for (let combo of sorted_combos) {
                    cumulativePct += combo.pct;
                    if (randomPick <= cumulativePct) {
                        selectedCombo = combo;
                        break;
                    }
                }
    
                if (selectedCombo) {
                    Decay.DecayInto(particle, selectedCombo.particles);
                }
            } else {
                console.warn(`No valid decay modes found for ${particle.name}`);
                return null;
            }
        }
        return null; // No decay
    }
    static DetermineCombinations(particle){
        let combos = [];
        let particles = Object.values(Constants.ParticleInfo);
        let num_particles = particles.length;
    
        for (let I = 0; I < num_particles; I++) {
            let pI = particles[I]; 
            for (let J = I + 1; J < num_particles; J++) {
                let pJ = particles[J];
                let m = pI.mass + pJ.mass;
                let c = pI.charge + pJ.charge;
                let b = pI.baryon + pJ.baryon;
                let le = pI.lepton.e + pJ.lepton.e;
                let lm = pI.lepton.m + pJ.lepton.m;
                let lt = pI.lepton.t + pJ.lepton.t;
    
                // 2-particle combinations
                let checkJ = Decay.CheckCombination(particle, m, c, b, le, lm, lt);
                if (checkJ === -1) break; // future pJ only increase mass
                if (checkJ === 1) combos.push({ particles: [Constants.ParticleInfo[pI.name], Constants.ParticleInfo[pJ.name]], pct: Math.pow((particle.mass - m), 2) });
    
                // 3-particle combinations
                for (let K = 0; K < num_particles; K++) {
                    let pK = particles[K];
                    let checkK = Decay.CheckCombination(particle,
                        m + pK.mass,
                        c + pK.charge,
                        b + pK.baryon,
                        le + pK.lepton.e,
                        lm + pK.lepton.m,
                        lt + pK.lepton.t
                    );
                    if (checkK === -1) break; // future pK only increase mass
                    if (checkK === 1) combos.push({ particles: [Constants.ParticleInfo[pI.name], Constants.ParticleInfo[pJ.name], Constants.ParticleInfo[pK.name]], pct: Math.pow((particle.mass - m), 2) });
    
                    // 4-particle combinations
                    for (let L = 0; L < num_particles; L++) {
                        let pL = particles[L];
                        let checkL = Decay.CheckCombination(particle,
                            m + pK.mass + pL.mass,
                            c + pK.charge + pL.charge,
                            b + pK.baryon + pL.baryon,
                            le + pK.lepton.e + pL.lepton.e,
                            lm + pK.lepton.m + pL.lepton.m,
                            lt + pK.lepton.t + pL.lepton.t
                        );
                        if (checkL === -1) break; // future pL only increase mass
                        if (checkL === 1) combos.push({ particles: [Constants.ParticleInfo[pI.name], Constants.ParticleInfo[pJ.name], Constants.ParticleInfo[pK.name], Constants.ParticleInfo[pL.name]], pct: Math.pow((particle.mass - m), 2) });
                    }
                }
            }
        }
    
        return combos;
    }
    
    static CheckCombination(p1, m, c, b, le, lm, lt){
        if(p1.mass > m){
            return (p1.charge === c &&
                    p1.baryon === b &&
                    p1.lepton.e === le &&
                    p1.lepton.m === lm &&
                    p1.lepton.t === lt) ? 1 : 0;
        }
        return -1;
    }
    
    static DecayInto(particle, decayProducts) {
        // 1. Check decay is possible (mass threshold)
        const totalDecayMass = decayProducts.reduce((sum, p) => sum + p.mass, 0);
        if (totalDecayMass > particle.mass) {
            console.error("Decay impossible: insufficient mass");
            return;
        }
    
        // 2. Calculate energies/momenta using conservation laws
        const E_parent = particle.mass; // Parent at rest
        let p = []; // Momenta magnitudes for decay products
        let directions = []; // Unit vectors for momentum directions
    
        // 3. Handle common decay types
        switch(decayProducts.length) {
            case 2: // Two-body decay (exact solution)
                const m1 = decayProducts[0].mass;
                const m2 = decayProducts[1].mass;
                
                const p_mag = Math.sqrt((E_parent**2 - (m1 + m2)**2) * (E_parent**2 - (m1 - m2)**2)) / (2 * E_parent);
    
                // Opposite directions for momentum conservation
                p = [p_mag, p_mag];
                directions = [
                    {x: 1, y: 0},  // Arbitrary direction
                    {x: -1, y: 0} // Opposite direction
                ];
                break;
    
            case 3: // Three-body decay (simplified phase space)
                // Use iterative method or precomputed tables
                // This is a placeholder - real implementation requires
                // Dalitz plot calculations or experimental data
                decayProducts.forEach(prod => {
                    const energy = (E_parent - totalDecayMass) * Math.random() + prod.mass;
                    p.push(Math.sqrt(energy**2 - prod.mass**2));
                });
                directions = Decay.GenerateIsotropicDirections(3); // Random angles
                break;
    
            default: // Generic multi-body (energy-momentum not strictly conserved)
                console.warn("Using simplified multi-body decay model");
                decayProducts.forEach(prod => {
                    const energy = E_parent/decayProducts.length;
                    p.push(Math.sqrt(energy**2 - prod.mass**2));
                });
                directions = Decay.GenerateIsotropicDirections(decayProducts.length);
        }
    
        // 4. Create particles with calculated momenta
        decayProducts.forEach((prod, i) => {
            const v = (p[i] / prod.mass) / Math.sqrt(1 + (p[i]/prod.mass)**2); // Relativistic speed
            const vx = v * directions[i].x;
            const vy = v * directions[i].y;
            Simulation.AddParticle(prod.name, particle.x, particle.y, vx, vy);
        });
    
        Simulation.RemoveParticle(particle);
    }
    
    // Helper for random direction vectors
    static GenerateIsotropicDirections(n) {
        return Array.from({length: n}, () => {
            const angle = 2 * Math.PI * Math.random();
            return {x: Math.cos(angle), y: Math.sin(angle)};
        });
    }

}

function applyBremsstrahlung({vx, vy, ax, ay, mass, charge, dt}) {
  const c = 3e8;
  const h = 6.626e-34;//plank constant

  const a2 = ax * ax + ay * ay;
  if (a2 === 0) return { new_vx: vx, new_vy: vy, photon: null };

  const epsilon0 = 8.854e-12;
  const P = (charge ** 2 * a2) / (6 * Math.PI * epsilon0 * c ** 3);
  const E_gamma = P * dt;

  // Photon direction (unit vector along acceleration)
  const a_mag = Math.sqrt(a2);
  const ax_unit = ax / a_mag;
  const ay_unit = ay / a_mag;

  // Photon momentum
  const p_gamma = E_gamma / c;
  const p_gamma_x = p_gamma * ax_unit;
  const p_gamma_y = p_gamma * ay_unit;

  // Updated particle momentum
  const p_x = mass * vx - p_gamma_x;
  const p_y = mass * vy - p_gamma_y;

  // New particle velocity
  const new_vx = p_x / mass;
  const new_vy = p_y / mass;

  // Photon wavelength
  const lambda = (h * c) / E_gamma;

  return {
    new_vx,
    new_vy,
    photon: {
      vx: ax_unit * c,
      vy: ay_unit * c,
      wavelength: lambda
    }
  };
}

