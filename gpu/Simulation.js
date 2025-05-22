class Simulation {
    static SIM = null;
    static PAUSED = true;
    static USER_INTERACTION_OCCURING = false;
    static SIM_STATS = {
        RUN_COUNT : 0,
        FRAME_COUNT : 0,
        FIRST_START : 0,
        LAST_RUN_START : 0,
        LAST_FRAME_START : 0,
        ELAPSED_RUN_TIME : 0,
        ELAPSED_FRAME_TIME : 0,
        ELAPSED_SIM_TIME : 0
    };
    static TIME_SCALAR = 1e-15;
    
    static async Init(particles){
        if(Simulation.SIM != null){
            throw new Error("Simulation Already Init'd");
        }
        if (!navigator.gpu) {
            throw new Error("WebGPU not supported");
        }
        
        const adapter = await navigator.gpu.requestAdapter();
        const device = await adapter.requestDevice();
        Simulation.SIM = new Simulation(device, particles);
    }

    static Run(after_frame_func = null){
        Simulation.AfterFrame = after_frame_func || ()=>{};
        if(Simulation.SIM){
            Simulation.SIM_STATS.FIRST_START = Date.now();
            requestAnimationFrame(Simulation._Run);
        }
    }
    
    static _Run(time_delta){
        Simulation.SIM_STATS.RUN_COUNT++;
        Simulation.SIM_STATS.ELAPSED_RUN_TIME+=time_delta;
        
        if(!Simulation.PAUSED && !Simulation.USER_INTERACTION_OCCURING){
            Simulation.SIM_STATS.FRAME_COUNT++;
            Simulation.SIM_STATS.ELAPSED_REAL_TIME += time_delta;
            time_delta *= Simulation.TIME_SCALAR;
            Simulation.SIM_STATS.ELAPSED_SIM_TIME += time_delta;
            
            Simulation.SIM.Update(time_delta);
            
            Simulation.AfterFrame();
        }
        requestAnimationFrame(Simulation._Run);
    }
    
    
    constructor(device, particles){
        this.device = device;
        this.particle_buffer = this.ToBuffer(particles);

        this.canvas = document.getElementById('canvas');
        this.context = this.canvas.getContext('webgpu');
        this.format = navigator.gpu.getPreferredCanvasFormat();
        this.context.configure({
            device: this.device,
            format: this.format,
            alphaMode: "opaque",
        });

        this.c = new Compute(device, this.particle_buffer);
        this.r = new Render(device, this.particle_buffer);
    }

    ToBuffer(particles){
        const particleBuffer = this.device.createBuffer({
            size: particleData.byteLength,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
            mappedAtCreation: true,
        });
        new Float32Array(particleBuffer.getMappedRange()).set(particleData);
        particleBuffer.unmap();
    }

    Update(delta_time){
        this.command_encoder = this.device.createCommandEncoder();
        this.UpdatePhysics();
        this.Render();
        this.device.queue.submit([this.command_encoder.finish()]);
    }

    UpdatePhysics(){
        const pass = commandEncoder.beginComputePass();
        pass.setPipeline(this.c.compute_pipeline);
        pass.setBindGroup(0, this.c.compute_bind_group);
        pass.dispatchWorkgroups(Math.ceil(PARTICLE_COUNT / 64));
        pass.end();
    }

    Render(){
        const textureView = context.getCurrentTexture().createView();
        const renderPass = this.command_encoder.beginRenderPass({
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
        renderPass.draw(6, PARTICLE_COUNT, 0, 0); // 6 verts per quad, N instances
        renderPass.end();
    }
}
