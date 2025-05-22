class Compute {
    constructor(device, particle_buffer, sim_uniform_buffer, particle_count) {
        this.device = device;
        this.particle_buffer = particle_buffer;
        this.sim_uniform_buffer = sim_uniform_buffer;
        this.particle_count = particle_count;

        this.compute_module = device.createShaderModule({ code: this.ComputeShader() });
        this.compute_pipeline = device.createComputePipeline({
            layout: 'auto',
            compute: { module: this.compute_module, entryPoint: 'main' },
        });

        this.particle_bind_group = device.createBindGroup({
            layout: this.compute_pipeline.getBindGroupLayout(0),
            entries: [{ binding: 0, resource: { buffer: this.particle_buffer } }],
        });
        this.sim_bind_group = device.createBindGroup({
            layout: this.compute_pipeline.getBindGroupLayout(1),
            entries: [{ binding: 0, resource: { buffer: this.sim_uniform_buffer } }],
        });
    }

    ComputeShader() {
        return `
            struct Particle {
                id: u32,x: f32,y: f32,vx: f32,vy: f32,mass: f32,charge: f32,symbol: u32,color: vec4<f32>,radius: f32,_pad: vec3<f32>,
            };
            
            struct SimOptions {
                dt : f32
            };
            
            @group(0) @binding(0) var<storage, read_write> particles: array<Particle>;
            @group(1) @binding(0) var<uniform> sim : SimOptions;
            
            const G: f32 = 6.67430e-11;
            const K: f32 = 8.9875517923e9;
            const C: f32 = 299792458.0;
            
            ${this.fn_compute_force()}
            ${this.fn_lorentz_gamma()}
            ${this.fn_velocity_from_momentum()}
            ${this.fn_clamp_velocity()}
            ${this.fn_particles_overlap()}
            ${this.fn_collision_response()}
            
            @compute @workgroup_size(64)
            fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
                let i = global_id.x;
                if (i >= arrayLength(&particles)) { return; }
            
                let pi = particles[i];
                var force = vec2<f32>(0.0, 0.0);
            
                // Sum all pairwise forces
                for (var j: u32 = 0u; j < arrayLength(&particles); j = j + 1u) {
                    if (i == j) { continue; }
                    let pj = particles[j];
                    force += compute_force(pi, pj);
                }
            
                // Relativistic momentum update
                let gamma = lorentz_gamma(pi.vx, pi.vy);
                var px = gamma * pi.mass * pi.vx + force.x * sim.dt;
                var py = gamma * pi.mass * pi.vy + force.y * sim.dt;
            
                // Get new velocity from momentum
                var v_new = velocity_from_momentum(px, py, pi.mass);
                v_new = clamp_velocity(v_new);

                var new_x = pi.x + v_new.x * sim.dt;
                var new_y = pi.y + v_new.y * sim.dt;
            
                // Collision handling (in-place)
                for (var j: u32 = 0u; j < arrayLength(&particles); j = j + 1u) {
                    if (i == j) { continue; }
                    let pj = particles[j];
                    v_new = collision_response(pi, pj, v_new, vec2<f32>(new_x, new_y));
                }
            
                // Update particle (symbol and color are preserved)
                particles[i].vx = v_new.x;
                particles[i].vy = v_new.y;
                particles[i].x = pi.x + v_new.x * sim.dt;
                particles[i].y = pi.y + v_new.y * sim.dt;
                // particles[i].symbol and particles[i].color remain unchanged
            }
        `;
    }

    fn_compute_force() {
        return `
            fn compute_force(pi: Particle, pj: Particle) -> vec2<f32> {
                let dx = pj.x - pi.x;
                let dy = pj.y - pi.y;
                let dist2 = dx * dx + dy * dy + 1e-6;
                let dist = sqrt(dist2);
            
                // Gravity
                let fg = G * pi.mass * pj.mass / dist2;
                // Electromagnetic
                let fem = K * pi.charge * pj.charge / dist2;
            
                // Total force vector
                let f_total = (fg + fem) / dist;
                return vec2<f32>(f_total * dx, f_total * dy);
            }`;
    }

    fn_lorentz_gamma() {
        return `
            fn lorentz_gamma(vx: f32, vy: f32) -> f32 {
                let v2 = vx * vx + vy * vy;
                let c2 = C * C;
                return 1.0 / sqrt(1.0 - v2 / c2);
            }`;
    }

    fn_velocity_from_momentum() {
        return `
            fn velocity_from_momentum(px: f32, py: f32, mass: f32) -> vec2<f32> {
                let c2 = C * C;
                let p2 = px * px + py * py;
                let m2c2 = mass * mass * c2;
                let denom = sqrt(m2c2 + p2);
                let vx = px * C / denom;
                let vy = py * C / denom;
                return vec2<f32>(vx, vy);
            }`;
    }

    fn_clamp_velocity() {
        return `
            fn clamp_velocity(v: vec2<f32>) -> vec2<f32> {
                let speed = length(v);
                if (speed > C) {
                    return v * (C / speed);
                }
                return v;
            }`;
    }

    fn_particles_overlap(){
        return `
            fn particles_overlap(x1: f32, y1: f32, x2: f32, y2: f32, min_dist: f32) -> bool {
                let dx = x2 - x1;
                let dy = y2 - y1;
                let dist2 = dx * dx + dy * dy;
                return dist2 < min_dist * min_dist;
            }`;
    }

    fn_collision_response(){
        return `
            fn collision_response(pi: Particle, pj: Particle,v_new: vec2<f32>, new_pos: vec2<f32>) -> vec2<f32> {
                let dx = pj.x - new_pos.x;
                let dy = pj.y - new_pos.y;
                let dist2 = dx * dx + dy * dy;
                let min_dist = pi.radius + pj.radius;
                if (dist2 < min_dist * min_dist) {
                    let dist = sqrt(dist2);
                    let nx = dx / dist;
                    let ny = dy / dist;
                    let rvx = v_new.x - pj.vx;
                    let rvy = v_new.y - pj.vy;
                    let vn = rvx * nx + rvy * ny;
                    if (vn < 0.0) {
                        let m1 = pi.mass;
                        let m2 = pj.mass;
                        let impulse = (2.0 * vn) / (m1 + m2);
                        let corr_vx = v_new.x - impulse * m2 * nx;
                        let corr_vy = v_new.y - impulse * m2 * ny;
                        // Optionally, push particle out of overlap (not updating position here)
                        return vec2<f32>(corr_vx, corr_vy);
                    }
                }
                return v_new;
            }`;
    }
}
