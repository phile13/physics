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
    ComputeShader(){
        return `
            struct Particle {
                id: u32, x: f32, y: f32, vx: f32, vy: f32, mass: f32, charge: f32, color: vec4<f32>,
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
            
                // Update particle
                particles[i].vx = v_new.x;
                particles[i].vy = v_new.y;
                particles[i].x = pi.x + v_new.x * sim.dt;
                particles[i].y = pi.y + v_new.y * sim.dt;
            }
        `;
    }

    fn_compute_force(){
        return `
            // Compute the force between two particles (gravity + EM)
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

    fn_lorentz_gamma(){
        return `
            // Compute Lorentz gamma factor
            fn lorentz_gamma(vx: f32, vy: f32) -> f32 {
                let v2 = vx * vx + vy * vy;
                let c2 = C * C;
                return 1.0 / sqrt(1.0 - v2 / c2);
            }`;
    }

    fn_velocity_from_momentum(){
        return `
            // Compute velocity from momentum (relativistic)
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

    fn_clamp_velocity(){
        return `
            // Clamp velocity to c if needed
            fn clamp_velocity(v: vec2<f32>) -> vec2<f32> {
                let speed = length(v);
                if (speed > C) {
                    return v * (C / speed);
                }
                return v;
            }`;
    }
    
}
