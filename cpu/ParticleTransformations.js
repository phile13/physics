class ParticleTransformations {

    static CheckForDecay(p, dt){
        const v_sqr = p.next.vx ** 2 + p.next.vy ** 2;
        const gamma = 1 / Math.sqrt(1 - Math.min(v_sqr / PhysicsEngine.C_SQR, 0.999999));
        p.age += dt / gamma;

        if(p.age >= p.life_span){
            return Decay.Run(p);
        }
        return null;
    }

    static CheckForFission(p, dt){
        if(Math.random() < p.fission_probability * dt){
            return Fission.Run(p);
        }
        return null;
    }

    static CheckForCollision(p1, p2, dt){
        if(p1.type === "Photon" && p2.type === "Photon" || p1.id === p2?.parent || p2.id === p1?.parent) return;

        const overlap = Particle.CheckForOverlap(p1,p2);
        if(overlap){
            if(p1.type === p2.anti_particle){
                const midx = 0.5 * (overlap.contactPoint1.x + overlap.contactPoint2.x);
                const midy = 0.5 * (overlap.contactPoint1.y + overlap.contactPoint2.y);
                // Example: produce two photons in random opposite directions
                const theta = Math.random() * 2 * Math.PI;
                const c = PhysicsEngine.C;
                const photon1 = new Particle("Photon", midx, midy, Math.cos(theta)*c, Math.sin(theta)*c);
                const photon2 = new Particle("Photon", midx, midy, -Math.cos(theta)*c, -Math.sin(theta)*c);
                return {
                    created: [photon1, photon2],
                    destroyed: [p1, p2]
                };
            }
            else if( false ){
                //hold for fusion
            }
            else if((p1.type === "Photon" && p2.type !== "Photon") || (p2.type === "Photon" && p1.type !== "Photon")){
                const photon = (p1.type === "Photon") ? p1 : p2;
                const target = (p1 === photon) ? p2 : p1;

                const c = PhysicsEngine.C;
                const pvx = photon.current.vx, pvy = photon.current.vy;
                const pSpeed = Math.hypot(pvx, pvy);
                const photonEnergy = photon.mass * c * c + photon.mass * pSpeed * pSpeed / 2;

                const rand = Math.random();

                // --- Absorption ---
                if (rand < 0.3) {
                    // Photon is absorbed, its energy is transferred to the target
                    target.next.vx += pvx * photon.mass / target.mass;
                    target.next.vy += pvy * photon.mass / target.mass;

                    return {
                        created: [],
                        destroyed: [photon]
                    };
                }

                // --- Compton-like Scattering ---
                else if (rand < 0.9) {
                    // Photon loses energy and changes direction
                    const scatterAngle = Math.random() * 2 * Math.PI;
                    const newSpeed = pSpeed * 0.8;

                    photon.next.vx = Math.cos(scatterAngle) * newSpeed;
                    photon.next.vy = Math.sin(scatterAngle) * newSpeed;

                    // Recoil: opposite direction to change in photon momentum
                    const deltaVx = (photon.next.vx - pvx);
                    const deltaVy = (photon.next.vy - pvy);

                    target.next.vx -= deltaVx * photon.mass / target.mass;
                    target.next.vy -= deltaVy * photon.mass / target.mass;

                    return null;
                }

                // --- Stimulated Emission ---
                else {
                    const emissionAngle = Math.random() * 2 * Math.PI;
                    const emittedPhoton = new Particle("Photon", target.current.x, target.current.y, Math.cos(emissionAngle) * c, Math.sin(emissionAngle) * c);
                    emittedPhoton.parent = target.id;

                    // Recoil to target
                    target.next.vx -= Math.cos(emissionAngle) * photon.mass / target.mass;
                    target.next.vy -= Math.sin(emissionAngle) * photon.mass / target.mass;

                    return {
                        created: [emittedPhoton],
                        destroyed: []
                    };
                }
            }
            else{//bounce
                const m1 = p1.mass, m2 = p2.mass;
                const v1x = p1.current.vx, v1y = p1.current.vy;
                const v2x = p2.current.vx, v2y = p2.current.vy;
                const nx = overlap.collisionNormal.x, ny = overlap.collisionNormal.y;
            
                // Relative velocity along the normal
                const relVel = (v1x - v2x) * nx + (v1y - v2y) * ny;
            
                // Don't resolve if not approaching
                if (relVel > 0) return null;
            
                // Impulse scalar (for 1D along collision normal)
                const impulse = (2 * relVel) / (m1 + m2);
            
                // Update velocities to conserve momentum along the normal
                p1.next.vx -= impulse * m2 * nx;
                p1.next.vy -= impulse * m2 * ny;
                p2.next.vx += impulse * m1 * nx;
                p2.next.vy += impulse * m1 * ny;
            }
        }
        return null;
    }

    static CheckForAnnihilation(p1, p2, dt){
        if(p1.type == p2.anti_particle){
            return null;
        }
        return null;
    }

    static CheckForFusion(p1, p2, dt){
        return null;
    }

    static CheckForHadronization(p1, p2, dt){
        return null;
    }

    static CheckForPairProduction(p1, p2, dt){
        return null;
    }

    static CheckForSynchrotronRadiation(p, dt) {
        if (p.is_force || p.charge === 0) return;
        const threshold = 1e-14; // minimum energy to emit photon
        const epsilon = 1e-12;

        // Momentum and gamma
        const p_sqr = p.px ** 2 + p.py ** 2;
        const gamma = Math.sqrt(1 + p_sqr / p.mass_sqr_times_c_sqr);
        const gamma_mass = gamma * p.mass;

        // Acceleration
        const ax = p.fx / gamma_mass;
        const ay = p.fy / gamma_mass;

        // Velocity
        const vx = p.current.vx;
        const vy = p.current.vy;

        // Perpendicular acceleration
        const v_dot_a = vx * ax + vy * ay;
        const v_mag_sqr = vx ** 2 + vy ** 2;

        const a_parallel_x = (v_dot_a / (v_mag_sqr + epsilon)) * vx;
        const a_parallel_y = (v_dot_a / (v_mag_sqr + epsilon)) * vy;

        const a_perp_x = ax - a_parallel_x;
        const a_perp_y = ay - a_parallel_y;
        const a_perp_sqr = a_perp_x ** 2 + a_perp_y ** 2;

        // Power
        const power = (p.charge ** 2 * a_perp_sqr * gamma ** 4) /
                      (6 * Math.PI * PhysicsEngine.EPSILON_0 * PhysicsEngine.C_CUBE);
        const energy_radiated = power * dt;

        if (energy_radiated > threshold) {
            const v_norm = Math.hypot(vx, vy) || 1;
            const dir_x = vx / v_norm;
            const dir_y = vy / v_norm;

            const photon = new Particle(
                "Photon",
                p.current.x,
                p.current.y,
                dir_x * PhysicsEngine.C,
                dir_y * PhysicsEngine.C
            );
            photon.energy = energy_radiated;
            photon.px = (energy_radiated / PhysicsEngine.C) * dir_x;
            photon.py = (energy_radiated / PhysicsEngine.C) * dir_y;
            photon.parent = p.id;
            return [photon];
        }
        return null;
    }
    

}

class Decay {
    static DECAY_TABLE = {
        Electron: null,
        Positron: null, // annihilation is a special case, not spontaneous decay
        Muon: ['Electron', 'Photon'],
        //Muon: ['Electron', 'ElectronAntineutrino', 'MuonNeutrino'],
        AntiMuon: ['Positron', 'MuonAntineutrino', 'ElectronNeutrino'],
        Tau: ['Electron', 'ElectronAntineutrino', 'TauNeutrino'], // or Muon version
        AntiTau: ['Positron', 'TauAntineutrino', 'ElectronNeutrino']
        // Add more...
    };

    static Run(p) {
        const products = Decay.DECAY_TABLE[p.type];
        if (!products) return null;

        const result = [];

        // Total momentum and energy of parent
        const total_px = p.px;
        const total_py = p.py;
        const total_E = Math.sqrt(p.px**2 + p.py**2 + p.mass_sqr_times_c_sqr) * PhysicsEngine.C;

        if (products.length === 2) {
            return this.GenerateTwoBodyDecay(products, total_px, total_py, total_E, p);
        }
        else if (products.length === 3){
            return this.GenerateThreeBodyDecay(products, total_px, total_py, total_E, p);
        }
        return result;
    }

    static PToV(px, py, m) {
        const c = PhysicsEngine.C;
        if (m === 0) {
            // Massless particle: velocity magnitude = c, direction = momentum direction
            const p_mag = Math.sqrt(px * px + py * py);
            if (p_mag === 0) return [0, 0]; // zero momentum edge case
            return [ (px / p_mag) * c, (py / p_mag) * c ];
        }
    
        // Massive particle: relativistic velocity
        const p_sqr = px * px + py * py;
        const gamma = Math.sqrt(1 + p_sqr / (m * m * c * c));
        const vx = px / (gamma * m);
        const vy = py / (gamma * m);
        return [vx, vy];
    }

    static GenerateTwoBodyDecay(types, parent_px, parent_py, parent_E, parent) {
        const [type1, type2] = types;
    
        const m1 = Particle.PARTICLE_TYPES[type1][1];
        const m2 = Particle.PARTICLE_TYPES[type2][1];
    
        const M = parent.mass;
        const E_cm = M * PhysicsEngine.C_SQR;
    
        const s = E_cm * E_cm;
        const m1c2 = m1 * PhysicsEngine.C_SQR;
        const m2c2 = m2 * PhysicsEngine.C_SQR;
    
        // Calculate magnitude of momentum in center-of-mass frame
        const p_mag = Math.sqrt(
            (s - (m1c2 + m2c2) ** 2) * (s - (m1c2 - m2c2) ** 2)
        ) / (2 * E_cm);
    
        // Random emission direction
        const angle = Math.random() * 2 * Math.PI;
        const dx = Math.cos(angle);
        const dy = Math.sin(angle);
    
        // Momentum vectors
        const px1 = p_mag * dx;
        const py1 = p_mag * dy;
        const px2 = -px1;
        const py2 = -py1;
    
        // Convert to velocities using relativistic formula
        const [vx1, vy1] = Decay.PToV(px1, py1, m1);
        const [vx2, vy2] = Decay.PToV(px2, py2, m2);
    
        // Position offset from center
        const r_offset = parent.radius || 1e-15;
        const x = parent.current.x;
        const y = parent.current.y;
    
        const x1 = x + dx * r_offset;
        const y1 = y + dy * r_offset;
        const x2 = x - dx * r_offset;
        const y2 = y - dy * r_offset;
    
        return [
            new Particle(type1, x1, y1, vx1, vy1),
            new Particle(type2, x2, y2, vx2, vy2)
        ];
    }
    

    static GenerateThreeBodyDecay(types, parent_px, parent_py, parent_E, parent) {
        const [t1, t2, t3] = types;
    
        const baseAngle = Math.random() * 2 * Math.PI;
        const angle1 = baseAngle;
        const angle2 = baseAngle + (2 * Math.PI / 3);
        const angle3 = baseAngle + (4 * Math.PI / 3);
    
        const pMag = 0.3 * parent.mass * PhysicsEngine.C;
    
        // Momentum vectors
        const px1 = pMag * Math.cos(angle1);
        const py1 = pMag * Math.sin(angle1);
        const px2 = pMag * Math.cos(angle2);
        const py2 = pMag * Math.sin(angle2);
        const px3 = -px1 - px2;
        const py3 = -py1 - py2;
    
        // Masses of decay products
        const m1 = Particle.PARTICLE_TYPES[t1][1];
        const m2 = Particle.PARTICLE_TYPES[t2][1];
        const m3 = Particle.PARTICLE_TYPES[t3][1];
    
        // Convert momentum to velocity
        const [vx1, vy1] = Decay.PToV(px1, py1, m1);
        const [vx2, vy2] = Decay.PToV(px2, py2, m2);
        const [vx3, vy3] = Decay.PToV(px3, py3, m3);
    
        // Radial position offset
        const r_offset = parent.radius || 1e-15;
        const x = parent.current.x;
        const y = parent.current.y;
    
        return [
            new Particle(t1, x + r_offset * Math.cos(angle1), y + r_offset * Math.sin(angle1), vx1, vy1),
            new Particle(t2, x + r_offset * Math.cos(angle2), y + r_offset * Math.sin(angle2), vx2, vy2),
            new Particle(t3, x + r_offset * Math.cos(angle3), y + r_offset * Math.sin(angle3), vx3, vy3)
        ];
    }
    
    
}
