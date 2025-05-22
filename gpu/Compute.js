class Compute {
    constructor(device, particle_buffer, particle_count){
        this.device = device;
        this.particle_buffer = particle_buffer;
        this.particle_count = particle_count;
        
        this.compute_module = device.createShaderModule({ code: this.ComputeShader() });
        this.compute_pipeline = device.createComputePipeline({
            layout: 'auto',
            compute: { module: this.compute_module, entryPoint: 'main' },
        });
        this.compute_bind_group = device.createBindGroup({
            layout: this.compute_pipeline.getBindGroupLayout(0),
            entries: [{ binding: 0, resource: { buffer: this.particle_buffer } }],
        });
    }

    ComputeShader(){
        return `
            struct Particle {
                id: u32, x: f32, y: f32, vx: f32, vy: f32, mass: f32, charge: f32,
            };
            @group(0) @binding(0) var<storage, read> particles: array<Particle>;
            
            @compute @workgroup_size(64)
            fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
                let i = global_id.x;
                if (i >= arrayLength(&particles)) { return; }
                
                // Simple gravity
                particles[i].vy = particles[i].vy + 0.1;
                
                // Update position
                particles[i].x = particles[i].x + particles[i].vx;
                particles[i].y = particles[i].y + particles[i].vy;
                
                // Boundary collision (canvas size 800x600)
                if (particles[i].x < 0.0) {
                    particles[i].x = 0.0;
                    particles[i].vx = -particles[i].vx * 0.8;
                }
                if (particles[i].x > 1000.0) {
                    particles[i].x = 1000.0;
                    particles[i].vx = -particles[i].vx * 0.8;
                }
                if (particles[i].y < 0.0) {
                    particles[i].y = 0.0;
                    particles[i].vy = -particles[i].vy * 0.8;
                }
                if (particles[i].y > 1000.0) {
                    particles[i].y = 1000.0;
                    particles[i].vy = -particles[i].vy * 0.8;
                }
            }
        `;
    }
}
