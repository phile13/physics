class Decay {
    static ParticleInfo = {
        Photon:        { name: "Photon",        lifespan: Infinity, charge: 0,  mass: 0,        baryon: 0, lepton: {e:0, m:0, t:0} },
        Gluon:         { name: "Gluon",         lifespan: Infinity, charge: 0,  mass: 0,        baryon: 0, lepton: {e:0, m:0, t:0} },
    
        NeutrinoE:     { name: "NeutrinoE",     lifespan: Infinity, charge: 0,  mass: 0.0000022,baryon: 0, lepton: {e:1, m:0, t:0} },
        AntiNeutrinoE: { name: "AntiNeutrinoE", lifespan: Infinity, charge: 0,  mass: 0.0000022,baryon: 0, lepton: {e:-1, m:0, t:0} },
 
        NeutrinoMu:    { name: "NeutrinoMu",    lifespan: Infinity, charge: 0,  mass: 0.17,     baryon: 0, lepton: {e:0, m:1, t:0} },
        AntiNeutrinoMu:{ name: "AntiNeutrinoMu",lifespan: Infinity, charge: 0,  mass: 0.17,     baryon: 0, lepton: {e:0, m:-1, t:0} },
    
        Electron:      { name: "Electron",      lifespan: Infinity, charge: -1, mass: 0.511,    baryon: 0, lepton: {e:1, m:0, t:0} },
        Positron:      { name: "Positron",      lifespan: Infinity, charge: +1, mass: 0.511,    baryon: 0, lepton: {e:-1, m:0, t:0} },
    
        Up:            { name: "Up",            lifespan: Infinity, charge: +2/3, mass: 2.2,    baryon: 1/3, lepton: {e:0, m:0, t:0} },
        AntiUp:        { name: "AntiUp",        lifespan: Infinity, charge: -2/3, mass: 2.2,    baryon: -1/3, lepton: {e:0, m:0, t:0} },
    
        Down:          { name: "Down",          lifespan: Infinity, charge: -1/3, mass: 4.7,    baryon: 1/3, lepton: {e:0, m:0, t:0} },
        AntiDown:      { name: "AntiDown",      lifespan: Infinity, charge: +1/3, mass: 4.7,    baryon: -1/3, lepton: {e:0, m:0, t:0} },
    
        NeutrinoTau:   { name: "NeutrinoTau",   lifespan: Infinity, charge: 0,  mass: 18.2,     baryon: 0, lepton: {e:0, m:0, t:1} },
        AntiNeutrinoTau:{name: "AntiNeutrinoTau",lifespan: Infinity, charge: 0, mass: 18.2,     baryon: 0, lepton: {e:0, m:0, t:-1} },
    
        Strange:       { name: "Strange",       lifespan: 1e-8,     charge: -1/3, mass: 96,     baryon: 1/3, lepton: {e:0, m:0, t:0} },
        AntiStrange:   { name: "AntiStrange",   lifespan: 1e-8,     charge: +1/3, mass: 96,     baryon: -1/3, lepton: {e:0, m:0, t:0} },
    
        Muon:          { name: "Muon",          lifespan: 2.2e-6,   charge: -1, mass: 105.7,    baryon: 0, lepton: {e:0, m:1, t:0} },
        AntiMuon:      { name: "AntiMuon",      lifespan: 2.2e-6,   charge: +1, mass: 105.7,    baryon: 0, lepton: {e:0, m:-1, t:0} },

        Pion:          { name: "Pion",          lifespan: 8.4e-17,   charge: 0, mass: 135.0,    baryon: 0, lepton: {e:0, m:0, t:0} },
        PionMinus:     { name: "PionMinus",     lifespan: 2.6e-8,   charge: -1, mass: 139.6,    baryon: 0, lepton: {e:0, m:0, t:0} },
        PionPlus:      { name: "PionPlus",      lifespan: 2.6e-8,   charge: +1, mass: 139.6,    baryon: 0, lepton: {e:0, m:0, t:0} },

        Charm:         { name: "Charm",         lifespan: 1e-12,    charge: +2/3, mass: 1280,   baryon: 1/3, lepton: {e:0, m:0, t:0} },
        AntiCharm:     { name: "AntiCharm",     lifespan: 1e-12,    charge: -2/3, mass: 1280,   baryon: -1/3, lepton: {e:0, m:0, t:0} },

        Tau:           { name: "Tau",           lifespan: 2.9e-13,  charge: -1, mass: 1776.9,   baryon: 0, lepton: {e:0, m:0, t:1} },
        AntiTau:       { name: "AntiTau",       lifespan: 2.9e-13,  charge: +1, mass: 1776.9,   baryon: 0, lepton: {e:0, m:0, t:-1} },
    
        Proton:        { name: "Proton",        lifespan: Infinity, charge: +1, mass: 938.3,    baryon: 1, lepton: {e:0, m:0, t:0} },
        AntiProton:    { name: "AntiProton",    lifespan: Infinity, charge: -1, mass: 938.3,    baryon: -1, lepton: {e:0, m:0, t:0} },

        Neutron:       { name: "Neutron",       lifespan: 880,      charge: 0,  mass: 939.6,    baryon: 1, lepton: {e:0, m:0, t:0} },
        AntiNeutron:   { name: "AntiNeutron",   lifespan: 880,      charge: 0,  mass: 939.6,    baryon: -1, lepton: {e:0, m:0, t:0} },

        Bottom:        { name: "Bottom",        lifespan: 1.5e-12,  charge: -1/3, mass: 4180,   baryon: 1/3, lepton: {e:0, m:0, t:0} },
        AntiBottom:    { name: "AntiBottom",    lifespan: 1.5e-12,  charge: +1/3, mass: 4180,   baryon: -1/3, lepton: {e:0, m:0, t:0} },
    
        WBosonPlus:    { name: "WBosonPlus",    lifespan: 3e-25, charge: +1, mass: 80379, baryon: 0, lepton: {e:0, m:0, t:0} },
        WBosonMinus:   { name: "WBosonMinus",   lifespan: 3e-25, charge: -1, mass: 80379, baryon: 0, lepton: {e:0, m:0, t:0} },
        ZBoson:        { name: "ZBoson",        lifespan: 3e-25,    charge: 0,  mass: 91187.6,  baryon: 0, lepton: {e:0, m:0, t:0} },
    
        Higgs:         { name: "Higgs",         lifespan: 1.6e-22,  charge: 0,  mass: 125100,   baryon: 0, lepton: {e:0, m:0, t:0} },

        Top:           { name: "Top",           lifespan: 5e-25,    charge: +2/3, mass: 173000, baryon: 1/3, lepton: {e:0, m:0, t:0} },
        AntiTop:       { name: "AntiTop",       lifespan: 5e-25,    charge: -2/3, mass: 173000, baryon: -1/3, lepton: {e:0, m:0, t:0} }
    };

    static DetermineCombinations(particle){
        let combos = [];
        let particles = Object.values(Decay.ParticleInfo);
        let num_particles = particles.length;

        for(let I = 0; I < num_particles; I++){
            let pI = particles[I]; 
            for(let J = I+1; J < num_particles; J++){
                let pJ = particles[J];
                let m = pI.mass + pJ.mass;
                let c = pI.charge + pJ.charge;
                let b = pI.baryon + pJ.baryon;
                let le = pI.lepton.e + pJ.lepton.e;
                let lm = pI.lepton.m + pJ.lepton.m;
                let lt = pI.lepton.t + pJ.lepton.t;

                let checkJ = Decay.CheckCombination(particle, m, c, b, le, lm, lt);
                if (checkJ === -1) break; // future pJ only increase mass
                if (checkJ === 1) combos.push({particles : [pI.name, pJ.name], pct : Math.pow((particle.mass - m), 2)});

                for(let K = 0; K < num_particles; K++){
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
                    if (checkK === 1) combos.push({particles : [pI.name, pJ.name, pK.name], pct : Math.pow((particle.mass - m), 2}));
                }
            }
        }
        return combos;
    }static HandleDecay(particle){
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
                // Create new particles based on the selected combination
                let decayProducts = selectedCombo.particles.map(name => {
                    let info = Decay.ParticleInfo[name];
                    return {
                        name: info.name,
                        charge: info.charge,
                        mass: info.mass,
                        baryon: info.baryon,
                        lepton: { ...info.lepton },
                        age: 0, // new particles are just born
                    };
                });

                return decayProducts;
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
    let particles = Object.values(Decay.ParticleInfo);
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

            let checkJ = Decay.CheckCombination(particle, m, c, b, le, lm, lt);
            if (checkJ === -1) break; // future pJ only increase mass
            if (checkJ === 1) combos.push({ particles: [Decay.ParticleInfo[pI.name], Decay.ParticleInfo[pJ.name]], pct: Math.pow((particle.mass - m), 2) });

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
                if (checkK === 1) combos.push({ particles: [Decay.ParticleInfo[pI.name], Decay.ParticleInfo[pJ.name], Decay.ParticleInfo[pK.name]], pct: Math.pow((particle.mass - m), 2) });

                // Now add the 4-particle combinations
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
                    if (checkL === 1) combos.push({ particles: [Decay.ParticleInfo[pI.name], Decay.ParticleInfo[pJ.name], Decay.ParticleInfo[pK.name], Decay.ParticleInfo[pL.name]], pct: Math.pow((particle.mass - m), 2) });
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
    
    static CheckForDecay(particle){
        return (particle.age > Decay.ParticleInfo[particle.name].lifespan);
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
                    Decay.CalculateDecay(particle, selectedCombo.particles);
                }
            } else {
                console.warn(`No valid decay modes found for ${particle.name}`);
                return null;
            }
        }
        return null; // No decay
    }

    // General decay calculation for 2 to 4 products
    static CalculateDecay(particle, decayProducts) {
        const totalEnergy = Decay.EnergyMomentum(particle);
        const totalMomentum = particle.momentum;

        // Determine the number of decay products (2, 3, or 4)
        const numDecayProducts = decayProducts.length;

        // Total mass of decay products
        let totalDecayMass = 0;
        decayProducts.forEach(prod => totalDecayMass += prod.mass);

        // Simplified approach: Calculate energies and momenta for the decay products
        // Note: More complex kinematics should be applied for realistic scenarios

        let decayResults = [];
        for (let i = 0; i < numDecayProducts; i++) {
            const mass = decayProducts[i].mass;
            const energy = totalEnergy / numDecayProducts; // Distribute energy equally (idealized)
            const momentum = Math.sqrt(Math.pow(energy, 2) - Math.pow(mass, 2));
            //todo add vx, vy
            Simulation.AddParticle(new Particle(decayProducts[i].name, particle.x, particle.y, vx, vy);
        }

        particle.destroy();
    }
}
