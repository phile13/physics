class PhysicsEngine {
    static G = -6.67430e-11;
    static K = 8.98755e9;
    static MU_0 = 4 * Math.PI * 1e-7;
    static MU_0_OVER_4PI = 1e-7; // Magnetic constant / 4π (μ₀/4π)
    static ONE_OVER_FOUR_PI_EPSILON_0 = (1 / (4 * Math.PI * 8.85418782e-12));
    static C = 3e8;
    static C_SQR = PhysicsEngine.C * PhysicsEngine.C;
    static C_CUBE = PhysicsEngine.C * PhysicsEngine.C * PhysicsEngine.C;
    static EPSILON_0 = 8.8541878128e-12;

    constructor(particle_event_callback) {
        this.particle_event_callback = particle_event_callback;
        this.time_manager = new TimeManager(10000);
        this.current_time = 0;
        this.frame = 0;
    }

    Init(particles){
        this.time_manager.Store(particles, 0);
        this.time_manager.Init(particles);
    }

    Update(particles, dt){
        this.ClearForces(particles);

        this.frame++;
        this.current_time = this.frame * dt;
        this.time_manager.Store(particles, this.current_time);


        this.ComputeForces(particles, dt);
        this.ApplyForces_SimplifiedVerlet(particles, dt);

        this.CheckDecays(particles, dt);
        this.CheckCollisions(particles, dt)
        // Optional: clear the cache if you want to avoid using old indices
        // this.clearRetardedTimeCache();
    }

    CleanupAfterDeletion(deletedParticles) {
        this.time_manager.Delete(deletedParticles);
    }



    ClearForces(particles){
        let num_particles = particles.length;
        for(let p = 0; p < num_particles; p++){
            particles[p].Clear();
        }
    }

    ComputeForces(particles, dt) {
        let num_particles = particles.length;
        for (let i = 0; i < num_particles; i++) {
            const p1 = particles[i];
            if (p1.is_force) continue;
    
            for (let j = 0; j < num_particles; j++) {
                if (i === j) continue;
                const p2 = particles[j];
                if (p2.is_force) continue;

                const offsets = this.time_manager.Get(p1, p2);
                if(offsets){
                    // Gravity
                    const Fg_mag = PhysicsEngine.G * p1.mass * p2.mass / offsets.r_sqr;
                    const Fg_x = Fg_mag * offsets.unit_x;
                    const Fg_y = Fg_mag * offsets.unit_y;
        
                    // Electric field (Coulomb, retarded)
                    const E_mag = PhysicsEngine.K * p2.charge / offsets.r_sqr;
                    const E_x = E_mag * offsets.unit_x;
                    const E_y = E_mag * offsets.unit_y;
        
                    // 
                    const Fe_x = p1.charge * E_x;
                    const Fe_y = p1.charge * E_y;
        
                    const B_z = (1.0 / PhysicsEngine.C_SQR) * (offsets.observed_data.vx * E_y - offsets.observed_data.vy * E_x);
                    const Fm_x = p1.charge * (p1.current.vy * B_z);
                    const Fm_y = p1.charge * (-p1.current.vx * B_z);


                    // Accumulate
                    p1.fx += Fg_x + Fe_x + Fm_x;
                    p1.fy += Fg_y + Fe_y + Fm_y;
                }
            }
        }
        
        this.CheckRadiation(particles, dt);
    }
    
    ApplyForces_SimplifiedVerlet(particles, dt) {
        for (const p of particles) {
            if (p.is_force) {
                // Copy velocity to next state, no force update
                p.next.vx = p.current.vx;
                p.next.vy = p.current.vy;
    
                // Update next position using next velocity
                p.next.x = p.current.x + p.next.vx * dt;
                p.next.y = p.current.y + p.next.vy * dt;
            } else {
                // Update relativistic momentum from force * dt
                p.px += p.fx * dt;
                p.py += p.fy * dt;
    
                // Calculate gamma * mass = sqrt(1 + p^2/(m^2 c^2)) * m
                const p_p_sqr = p.px * p.px + p.py * p.py;
                const gamma_times_mass = Math.sqrt(1 + p_p_sqr / p.mass_sqr_times_c_sqr) * p.mass;
    
                // Update next velocity from momentum
                p.next.vx = p.px / gamma_times_mass;
                p.next.vy = p.py / gamma_times_mass;
    
                // Update next position using current next velocity + approx acceleration term
                // Approximate acceleration = force / (gamma * mass)
                const ax = p.fx / gamma_times_mass;
                const ay = p.fy / gamma_times_mass;
    
                p.next.x = p.current.x + p.next.vx * dt + 0.5 * ax * dt * dt;
                p.next.y = p.current.y + p.next.vy * dt + 0.5 * ay * dt * dt;
            }
        }
    }

    CheckRadiation(particles, dt){
        for (let i = particles.length - 1; i >= 0; i--) {
            const p = particles[i];
            if (p.is_force) continue;

            const radiated = ParticleTransformations.CheckForSynchrotronRadiation(p, dt);
            if (radiated) {
                if(this.frame % 10 == 0){
                    this.particle_event_callback([], radiated);
                }
            }
        }
    }

    CheckDecays(particles, dt){
        for (let i = particles.length - 1; i >= 0; i--) {
            const p = particles[i];
            if (p.is_force) continue;

            const decayed = ParticleTransformations.CheckForDecay(p, dt);
            if (decayed) {
                this.particle_event_callback([p], decayed);
            }

            const fissioned = ParticleTransformations.CheckForFission(p, dt);
            if(fissioned){
                this.particle_event_callback([p], fissioned);
            }
        }
    }

    CheckCollisions(particles, dt){
        let num_particles = particles.length;
        for (let i = 0; i < num_particles; i++) {
            const p1 = particles[i];    
            for (let j = i+1; j < num_particles; j++) {
                const p2 = particles[j];
                
                const collision = ParticleTransformations.CheckForCollision(p1, p2);
                if(collision){
                    this.particle_event_callback(collision.destroyed, collision.created);
                }
            }
        }
    }

    CheckCompsite(particles, dt){
        // TODO hadronization
    }

    CheckPhotons(){
        // TODO pair production
    }

    DrawTrajectories(ctx, transformer) {
        this.time_manager.DrawPositionHistories(ctx, transformer);
    }

    DrawForcesHistory(ctx, transformer) {
        this.time_manager.DrawForcesHistory(ctx, transformer);
    }

    static GetPositionalDifferences(p1, p2){
        const dx = p1.current.x - p2.current.x;
        const dy = p1.current.y - p2.current.y;
        const r_sqr = Math.max(dx * dx + dy * dy , 1e-100);
        const r = Math.sqrt(r_sqr);
        return [dx,dy,r_sqr,r,dx/r,dy/r];
    }

    static GetVelocity(p){
        const p_sqr  = p.px ** 2 + p.py ** 2;
        const gamma_times_mass = Math.sqrt(1 + (p_sqr / p.mass_sqr_times_c_sqr)) * p.mass;
        const vx = p.px / gamma_times_mass;
        const vy = p.py / gamma_times_mass;
        return [vx, vy];
    }

    static GetAcceleration(p) {
        if (!p.prev_px || !p.prev_py) return [0, 0];
    
        const ax = (p.px - p.prev_px) / p.last_dt;
        const ay = (p.py - p.prev_py) / p.last_dt;
    
        return [ax, ay];
    }

    static AddBrakingForce(p, dt){
        const q = particle.charge;
        if(q === 0) return;

        const c = PhysicsEngine.C;
        const mu0 = 4 * Math.PI * 1e-7;

        const [vx, vy] = PhysicsEngine.getVelocity(particle);
        const vMag = Math.sqrt(vx * vx + vy * vy);
        if (vMag < 1e-8) return;

        // Compute acceleration using finite difference
        const [ax, ay] = PhysicsEngine.getAcceleration(particle);
        const aMag2 = ax * ax + ay * ay;

        // Radiated power (non-relativistic approximation)
        const power = (mu0 * q * q * aMag2) / (6 * Math.PI * PhysicsEngine.C);

        // Energy loss = power * dt
        const dE = power * dt;

        // Convert to momentum loss along velocity vector
        const pLoss = dE / PhysicsEngine.C;
        const fx = -pLoss * (vx / vMag) / dt;
        const fy = -pLoss * (vy / vMag) / dt;

        // Apply force
        particle.fx += fx;
        particle.fy += fy;
        particle.radiatedEnergy = (particle.radiatedEnergy || 0) + dE;
    }
}
