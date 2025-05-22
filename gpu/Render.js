class Render {

    constructor(device, particle_buffer){
        this.device = device;
        this.particle_buffer = particle_buffer;
        this.vertex_module =  device.createShaderModule({ code: this.VertexShader() });
        this.fragment_module = device.createShaderModule({ code: this.FragmentShader() });

        this.quad_verts = new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]);
        this.quad_buffer = device.createBuffer({
            size: this.quad_verts.byteLength,
            usage: GPUBufferUsage.VERTEX,
            mappedAtCreation: true,
        });
        new Float32Array(this.quad_buffer.getMappedRange()).set(this.quad_verts);
        this.quad_buffer.unmap();

        this.render_pipeline = device.createRenderPipeline({
            layout: 'auto',
            vertex: {
                module: this.vertex_module,
                entryPoint: 'main',
                buffers: [{
                    arrayStride: 2 * 4,
                    attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x2' }],
                }],
            },
            fragment: {
                module: this.fragment_module,
                entryPoint: 'main',
                targets: [{ format }],
            },
            primitive: { topology: 'triangle-list' },
        });

        const this.render_bind_group = device.createBindGroup({
            layout: this.render_pipeline.getBindGroupLayout(0),
            entries: [{ binding: 0, resource: { buffer: this.particle_buffer } }],
        });
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

    FragmentShader(){
        return `
            @fragment
            fn main(@location(0) local: vec2<f32>) -> @location(0) vec4<f32> {
                if (length(local) > 1.0) { discard; }
                return vec4<f32>(1.0, 1.0, 1.0, 1.0);
            }`;
    }

}
