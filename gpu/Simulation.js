class Simulation {
    static async Init(particles){
        if (!navigator.gpu) {
            throw new Error("WebGPU not supported");
        }
        const adapter = await navigator.gpu.requestAdapter();
        const device = await adapter.requestDevice();
        
        
    }
    
    constructor(){
        this.device = device;
        this.c = c;
        this.r = r;
        
        this.canvas = document.getElementById('canvas');
        this.context = this.canvas.getContext('webgpu');
        this.format = navigator.gpu.getPreferredCanvasFormat();
        context.configure({
            this.device,
            this.format,
            alphaMode: "opaque",
        });
    }

    Update(){
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
