class Renderer {
    constructor(device, particle_buffer, render_uniform_buffer, format, particle_count, atlasTexture, atlasSampler, atlasCols) {
        this.device = device;
        this.particle_buffer = particle_buffer;
        this.render_uniform_buffer = render_uniform_buffer;
        this.format = format;
        this.particle_count = particle_count;
        this.atlasTexture = atlasTexture;
        this.atlasSampler = atlasSampler;
        this.atlasCols = atlasCols; // number of columns in atlas grid

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
                { binding: 2, resource: this.atlasTexture.createView() },
                { binding: 3, resource: this.atlasSampler },
            ],
        });
    }

    VertexShader() {
        return `
            struct Particle {
                id: u32,x: f32,y: f32,vx: f32,vy: f32,mass: f32,charge: f32,symbol: u32,color: vec4<f32>,
            };
            struct RenderOptions {
                ulx: f32, uly: f32, inc: f32, pxw: u32, pxh: u32, showSymbol: u32, _pad1: u32,
            };
            
            @group(0) @binding(0) var<storage, read> particles: array<Particle>;
            @group(0) @binding(1) var<uniform> renderOpts : RenderOptions;
            
            struct VertexOut {
                @builtin(position) pos: vec4<f32>,
                @location(0) color: vec4<f32>,
                @location(1) local: vec2<f32>,
                @location(2) symbol: u32,
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
                let r_px = 16.0;
                let offset_px = localPos * r_px;
            
                // Offset in NDC
                let offset_ndc_x = offset_px.x / f32(renderOpts.pxw) * 2.0;
                let offset_ndc_y = -offset_px.y / f32(renderOpts.pxh) * 2.0;
            
                var out: VertexOut;
                out.pos = vec4<f32>(ndc_x + offset_ndc_x, ndc_y + offset_ndc_y, 0.0, 1.0);
                out.color = p.color;
                out.local = (localPos + vec2<f32>(1.0, 1.0)) * 0.5; // Map from [-1,1] to [0,1] for UV
                out.symbol = p.symbol;
                return out;
            }`;
    }

    FragmentShader() {
        // Atlas columns must be passed as a constant or uniform (here hardcoded for simplicity)
        return `
            @group(0) @binding(2) var atlas: texture_2d<f32>;
            @group(0) @binding(3) var atlasSampler: sampler;
            
            const ATLAS_COLS: u32 = 8u; // Set this to your atlas grid columns
            
            @fragment
            fn main(@location(0) color: vec4<f32>,@location(1) local: vec2<f32>,@location(2) symbol: u32) -> @location(0) vec4<f32> {
            
                if (length(local * 2.0 - 1.0) > 1.0) { discard; }

                if (renderOpts.showSymbol == 1u) {
                    let cellSize = 1.0 / f32(ATLAS_COLS);
                    let col = symbol % ATLAS_COLS;
                    let row = symbol / ATLAS_COLS;
                    let uv = vec2<f32>(
                        (f32(col) + local.x) * cellSize,
                        (f32(row) + local.y) * cellSize
                    );
                    let glyph = textureSample(atlas, atlasSampler, uv);
                    // Multiply glyph alpha by particle color
                    return vec4<f32>(color.rgb, color.a * glyph.a);
                }
                else{
                    return color;
                }
            }`;
    }

    render(passEncoder) {
        passEncoder.setPipeline(this.render_pipeline);
        passEncoder.setBindGroup(0, this.render_bind_group);
        passEncoder.setVertexBuffer(0, this.quad_buffer);
        passEncoder.draw(6, this.particle_count, 0, 0);
    }
}
