const { defineComponent } = Vue;


export const PixiPhysicsViewer = defineComponent({
    name: 'PixiPhysicsViewer',
    data() {
        return {
            container : null,
            canvas : null,
            isAnimating : true,
            fps: 0,
            frameCount: 0,
            simulatedTime : 0,
            timeScale : 1.0,
            xDimensions : 0,
            yDimensions : 0
        };
    },
    methods: {
        toggleAnimation() {
            this.isAnimating = !this.isAnimating;
            Simulation.UpdateTimeScale(this.isAnimating ? 1.0 : 0.0);
        },
        zoomToAll() {
            Simulation.AUTO_ZOOM_ON = !Simulation.AUTO_ZOOM_ON;
            //Simulation.AddToZoomParticle("bunny"); // or all particle IDs
        },
        addRandomParticle() {
            const id = `p${Date.now()}`;
            const radius = 5e-17;
            const type = 'e-';
            const x = Math.random() * 1e-15;
            const y = Math.random() * 1e-15;
            const vx = 0;
            const vy = 0;

            const p = new CircleParticle(id, type, x, y, vx, vy, radius);
            Simulation.AddParticle(p);
        },
        updateTimeScale(e) {
            this.timeScale = Math.pow(10, ((+e.target.value)/100) * (Math.log10(1e20) - Math.log10(1e-30)) + Math.log10(1e-30));
            Simulation.UpdateTimeScale(this.timeScale); // Use the time scale value from the slider
        }
    },
    computed: {
        formattedTimeScale() {
            return (this.timeScale) ? (+this.timeScale).toExponential(2) : 0; // Format timeScale to scientific notation with 2 decimals
        }
    },
    mounted() {
        Simulation.LinkPixi('simulation-canvas');
        this._fpsInterval = setInterval(() => {
            this.fps = Simulation.FRAME_STATS.fps.toFixed(2);
            this.frameCount = Simulation.FRAME_STATS.count;
            this.simulatedTime  = Simulation.FRAME_STATS.simulated_time.toFixed(2);
            this.xDimensions = (Simulation.SCREEN_INFO.WIDTH).toExponential(2);
            this.yDimensions = (Simulation.SCREEN_INFO.HEIGHT).toExponential(2);
        }, 500);
    },
    beforeUnmount() {
        clearInterval(this._fpsInterval);
    },
    template: `
        <div>
            <div ref="container" class="simulation-container">
                <canvas ref="canvas" id="simulation-canvas"></canvas>

                <div class="controls">
                    <button @click="toggleAnimation">{{ isAnimating ? 'Pause' : 'Play' }}</button>
                    <button @click="zoomToAll">Auto Zoom</button>
                    <button @click="addRandomParticle">Add Particle</button>
                </div>
                <div>
                    <label for="timeScaleSlider">Time Scale: {{ formattedTimeScale }}</label>
                    <input 
                        type="range" 
                        id="timeScaleSlider" 
                        min="0" 
                        max="100" 
                        step=".01"
                        @input="updateTimeScale" />
                </div>
                <div class="stats">
                    <p>FPS: {{ fps }}</p>
                    <p>Frame Count: {{ frameCount }}</p>
                    <p>Sim Time: {{ simulatedTime  }}</p>
                </div>
                <div>
                    X : {{ xDimensions }} -- Y : {{ yDimensions }} 
                </div>
            </div>
        </div>
    `
});
