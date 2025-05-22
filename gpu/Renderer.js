class Renderer {
    constructor(device, particle_buffer, render_uniform_buffer, format, particle_count) {
        this.device = device;
        this.particle_buffer = particle_buffer;
        this.render_uniform_buffer = render_uniform_buffer;
        this.format = format;
        this.particle_count = particle_count;

        // Full-screen quad for each particle (6 verts per instance)
        this.quad_verts = new Float32Array([
            -1, -1,  1, -1,  -1,  1,
            -1,  1,  1, -1,   1,  1
        ]);
        this.quad_buffer = device.createBuffer({
            size: this.quad_verts.byteLength,
            usage: GPUBufferUsage.VERTEX,
            mappedAtCreation: true,
        });
        new Float32Array(this.quad_buffer.getMappedRange()).set(this.quad_verts);
        this.quad_buffer.unmap();

        this.vertex_module = device.createShaderModule({ code: this.VertexShader() });
        this.fragment_module = device.createShaderModule({ code: this.FragmentShader() });

        this.render_pipeline = device.createRenderPipeline({
            layout: 'auto',
            vertex: {
                module: this.vertex_module,
                entryPoint: 'main',
                buffers: [
                    {
                        // Quad verts
                        arrayStride: 2 * 4,
                        attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x2' }],
                    }
                ],
            },
            fragment: {
                module: this.fragment_module,
                entryPoint: 'main',
                targets: [{ format: this.format }],
            },
            primitive: { topology: 'triangle-list' },
        });

        this.render_bind_group = device.createBindGroup({
            layout: this.render_pipeline.getBindGroupLayout(0),
            entries: [
                { binding: 0, resource: { buffer: this.particle_buffer } },
                { binding: 1, resource: { buffer: this.render_uniform_buffer } },
            ],
        });
    }

    VertexShader() {
        return `
            struct Particle {
                id: u32, x: f32, y: f32, vx: f32, vy: f32, mass: f32, charge: f32,
            };
            struct RenderOptions {
                ulx: f32, uly: f32, inc: f32, pxw: u32, pxh: u32, _pad0: u32, _pad1: u32,
            };
            
            @group(0) @binding(0) var<storage, read> particles: array<Particle>;
            @group(0) @binding(1) var<uniform> renderOpts : RenderOptions;
            
            struct VertexOut {
                @builtin(position) pos: vec4<f32>,
                @location(0) local: vec2<f32>,
            };
            
            @vertex
            fn main(@builtin(instance_index) instance: u32, @location(0) localPos: vec2<f32> ) -> VertexOut {
                let p = particles[instance];
            
                // Convert world to pixel
                let px = (p.x - renderOpts.ulx) / renderOpts.inc;
                let py = (p.y - renderOpts.uly) / renderOpts.inc;
            
                // Convert pixel to NDC
                let ndc_x = (px / f32(renderOpts.pxw)) * 2.0 - 1.0;
                let ndc_y = 1.0 - (py / f32(renderOpts.pxh)) * 2.0;
            
                // Particle radius in pixels
                let r_px = 4.0;
                let offset_px = localPos * r_px;
            
                // Offset in NDC
                let offset_ndc_x = offset_px.x / f32(renderOpts.pxw) * 2.0;
                let offset_ndc_y = -offset_px.y / f32(renderOpts.pxh) * 2.0;
            
                var out: VertexOut;
                out.pos = vec4<f32>(ndc_x + offset_ndc_x, ndc_y + offset_ndc_y, 0.0, 1.0);
                out.local = localPos;
                return out;
            }`;
    }

    FragmentShader() {
        return `
            @fragment
            fn main(@location(0) local: vec2<f32>) -> @location(0) vec4<f32> {
                if (length(local) > 1.0) { discard; }
                return vec4<f32>(1.0, 1.0, 1.0, 1.0);
            }`;
    }

    render(passEncoder) {
        passEncoder.setPipeline(this.render_pipeline);
        passEncoder.setBindGroup(0, this.render_bind_group);
        passEncoder.setVertexBuffer(0, this.quad_buffer);
        passEncoder.draw(6, this.particle_count, 0, 0); // 6 verts per quad, N instances
    }
}
