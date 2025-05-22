import { Compute } from './Compute.js';
import { Renderer } from './Renderer.js';

class Simulation {
    static SIM = null;
    static PAUSED = true;
    static USER_INTERACTION_OCCURING = false;
    static SHOW_SYMBOLS = false;
    static AUTO_ZOOM = false;
    static SIM_STATS = {
        RUN_COUNT: 0,
        FRAME_COUNT: 0,
        FIRST_START: 0,
        LAST_RUN_START: 0,
        LAST_FRAME_START: 0,
        ELAPSED_RUN_TIME: 0,
        ELAPSED_FRAME_TIME: 0,
        ELAPSED_SIM_TIME: 0,
        ELAPSED_REAL_TIME: 0
    };
    static TIME_SCALAR = 1e-15;

    // List of all particle/antiparticle symbols to be rendered in the atlas
    static SYMBOLS = [
        "e", "μ", "τ", "ν", "u", "d", "c", "s", "t", "b", "γ", "g", "W", "Z", "H", "p", "n",
        "Λ", "Σ", "Ξ", "Ω", "π",
        "e⁺", "e⁻", "μ⁺", "μ⁻", "π⁺", "π⁰", "π⁻", "K⁺", "K⁰", "K⁻", "Σ⁺", "Σ⁰", "Σ⁻", "Ξ⁰", "Ξ⁻", "Ω⁻",
        "p̄", "n̄", "Λ̄"
    ];

    static async Init(particles, symbols = null) {
        if (Simulation.SIM != null) {
            throw new Error("Simulation Already Init'd");
        }
        if (!navigator.gpu) {
            throw new Error("WebGPU not supported");
        }

        const adapter = await navigator.gpu.requestAdapter();
        const device = await adapter.requestDevice();

        // 1. Create the symbol atlas texture
        const { atlasTexture, symbolMap } = await Simulation.CreateSymbolAtlas(device, symbols || Simulation.SYMBOLS);

        // 2. Create uniform buffers
        const render_uniform_buffer = device.createBuffer({
            size: 32,
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
        });
        const sim_uniform_buffer = device.createBuffer({
            size: 4,
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
        });
        const energy_buffer = device.createBuffer({
            size: 4,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC,
        });
        const staging_buffer = device.createBuffer({
            size: 4,
            usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ,
        });

        // 3. Add symbol index to each particle
        particles.forEach(p => {
            // p.symbol should be the index in SYMBOLS array
            p.symbol = symbolMap[p.symbolName]; // e.g., p.symbolName = "π⁺"
        });

        Simulation.SIM = new Simulation(
            device,
            particles,
            render_uniform_buffer,
            sim_uniform_buffer,
            atlasTexture
        );
    }

    static async CreateSymbolAtlas(device, symbols) {
        // Create a canvas and draw each symbol in a grid
        const cellSize = 48;
        const cols = 8;
        const rows = Math.ceil(symbols.length / cols);
        const canvas = document.createElement('canvas');
        canvas.width = cols * cellSize;
        canvas.height = rows * cellSize;
        const ctx = canvas.getContext('2d');
        ctx.font = "bold 36px 'STIX Two', 'DejaVu Sans', 'Arial Unicode MS', sans-serif";
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";

        // Map symbol string to atlas index
        const symbolMap = {};
        symbols.forEach((sym, i) => {
            const col = i % cols;
            const row = Math.floor(i / cols);
            ctx.clearRect(col * cellSize, row * cellSize, cellSize, cellSize);
            ctx.fillText(sym, col * cellSize + cellSize / 2, row * cellSize + cellSize / 2);
            symbolMap[sym] = i;
        });

        // Upload canvas to GPU as a texture
        const imageBitmap = await createImageBitmap(canvas);
        const atlasTexture = device.createTexture({
            size: [canvas.width, canvas.height, 1],
            format: "rgba8unorm",
            usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT,
        });
        device.queue.copyExternalImageToTexture(
            { source: imageBitmap },
            { texture: atlasTexture },
            [canvas.width, canvas.height]
        );

        return { atlasTexture, symbolMap };
    }

    static Run(after_frame_func = null) {
        Simulation.AfterFrame = after_frame_func || ((x) => {});
        if (Simulation.SIM) {
            Simulation.SIM_STATS.FIRST_START = Date.now();
            requestAnimationFrame(Simulation._Run);
        }
    }

    static _Run(time_delta) {
        Simulation.SIM_STATS.RUN_COUNT++;
        Simulation.SIM_STATS.ELAPSED_RUN_TIME += time_delta;

        if (!Simulation.PAUSED && !Simulation.USER_INTERACTION_OCCURING) {
            Simulation.SIM_STATS.FRAME_COUNT++;
            Simulation.SIM_STATS.ELAPSED_REAL_TIME += time_delta;
            const scaled_dt = time_delta * Simulation.TIME_SCALAR;
            Simulation.SIM_STATS.ELAPSED_SIM_TIME += scaled_dt;

            // Update simulation uniform (dt)
            Simulation.SIM.device.queue.writeBuffer(
                Simulation.SIM.sim_uniform_buffer,
                0,
                new Float32Array([scaled_dt])
            );

            Simulation.SIM.Update(scaled_dt);
        }
        Simulation.AfterFrame(Simulation.SIM_STATS);
        requestAnimationFrame(Simulation._Run);
    }

    constructor(device, particles, render_uniform_buffer, sim_uniform_buffer, atlasTexture) {
        this.device = device;
        this.particle_count = particles.length;
        this.particle_buffer = this.ToBuffer(particles);
        this.render_uniform_buffer = render_uniform_buffer;
        this.sim_uniform_buffer = sim_uniform_buffer;
        this.atlasTexture = atlasTexture;

        this.canvas = document.getElementById('canvas');
        this.context = this.canvas.getContext('webgpu');
        this.format = navigator.gpu.getPreferredCanvasFormat();
        this.context.configure({
            device: this.device,
            format: this.format,
            alphaMode: "opaque",
        });

        this.c = new Compute(
            device,
            this.particle_buffer,
            this.sim_uniform_buffer,
            this.particle_count
        );
        this.r = new Render(
            device,
            this.particle_buffer,
            this.render_uniform_buffer,
            this.format,
            this.particle_count,
            this.atlasTexture,
            Simulation.SYMBOLS.length // pass atlas grid info as needed
        );
    }

    ToBuffer(particles) {
        // Each particle: id, x, y, vx, vy, mass, charge, symbol, r, g, b, a
        const particleData = new Float32Array(particles.flatMap(p => [
            p.id, p.x, p.y, p.vx, p.vy, p.mass, p.charge, p.symbol,
            p.color.r, p.color.g, p.color.b, p.color.a
        ]));

        const particleBuffer = this.device.createBuffer({
            size: particleData.byteLength,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
            mappedAtCreation: true,
        });
        new Float32Array(particleBuffer.getMappedRange()).set(particleData);
        particleBuffer.unmap();
        return particleBuffer;
    }

    Update(delta_time) {
        const commandEncoder = this.device.createCommandEncoder();

        // Compute pass
        const computePass = commandEncoder.beginComputePass();
        computePass.setPipeline(this.c.compute_pipeline);
        computePass.setBindGroup(0, this.c.particle_bind_group);
        computePass.setBindGroup(1, this.c.sim_bind_group);
        computePass.dispatchWorkgroups(Math.ceil(this.particle_count / 64));
        computePass.end();

        // Render pass
        const textureView = this.context.getCurrentTexture().createView();
        const renderPass = commandEncoder.beginRenderPass({
            colorAttachments: [{
                view: textureView,
                clearValue: { r: 0.07, g: 0.07, b: 0.07, a: 1 },
                loadOp: 'clear',
                storeOp: 'store',
            }],
        });
        renderPass.setPipeline(this.r.render_pipeline);
        renderPass.setBindGroup(0, this.r.render_bind_group);
        renderPass.setVertexBuffer(0, this.r.quad_buffer);
        renderPass.draw(6, this.particle_count, 0, 0);
        renderPass.end();

        this.device.queue.submit([commandEncoder.finish()]);
    }

    function updateRenderUniformBuffer({ulx, uly, inc, pxw, pxh, showSymbol}) {
        const buffer = new ArrayBuffer(32);
        const f32 = new Float32Array(buffer);
        const u32 = new Uint32Array(buffer);
        f32[0] = ulx;
        f32[1] = uly;
        f32[2] = inc;
        u32[3] = pxw;
        u32[4] = pxh;
        u32[5] = showSymbol ? 1 : 0; // NEW: showSymbol
        u32[6] = 0; // padding
        device.queue.writeBuffer(renderUniformBuffer, 0, buffer);
    }

}
