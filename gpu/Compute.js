class Compute {
    constructor(){

    }

    ComputeShader(){
        return `
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
