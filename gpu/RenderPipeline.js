class RenderPipeline {

    constructor(device){
        this.device = device;
        this.vertex_module = 
    }

    VertexShader(){
        return `
            struct Particle {
                id: u32, x: f32, y: f32, vx: f32, vy: f32, mass: f32, charge: f32,
            };
            @group(0) @binding(0) var<storage, read> particles: array<Particle>;
            
            struct VertexOut {
                @builtin(position) pos: vec4<f32>,
                @location(0) local: vec2<f32>,
            };
            
            @vertex
            fn main(@builtin(instance_index) instance: u32, @location(0) localPos: vec2<f32>) -> VertexOut {
                let p = particles[instance];
                // Circle radius
                let r = 4.0;
                let world = vec2<f32>(p.x, p.y) + localPos * r;
                // Convert to NDC
                let ndc = vec2<f32>(world.x / 500.0 - 1.0, 1.0 - world.y / 500.0);
                var out: VertexOut;
                out.pos = vec4<f32>(ndc, 0.0, 1.0);
                out.local = localPos;
                return out;
            }`;
    }

}
