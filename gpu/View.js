// View.js
// Handles pan/zoom, HUD, user interaction, and toggle button events.

export class View {
    constructor(canvas, simulation) {
        this.canvas = canvas;
        this.simulation = simulation;

        // Pan/zoom state
        this.pan = { x: -canvas.width / 2, y: -canvas.height / 2 };
        this.zoom = 1;
        this.minZoom = 0.1;
        this.maxZoom = 10;

        // Mouse interaction state
        this.isDragging = false;
        this.lastMouse = { x: 0, y: 0 };
        this.wheelDebounceTimer = null;

        // HUD overlay
        this._createHud();

        // Bind user interaction and control events
        this._bindEvents();
        this._bindControlEvents();

        // Set initial render options
        this.updateRenderOptions();
    }

    _createHud() {
        this.hud = document.createElement('div');
        Object.assign(this.hud.style, {
            position: 'absolute',
            left: '10px',
            top: '10px',
            color: '#0f0',
            background: 'rgba(0,0,0,0.7)',
            font: '16px monospace',
            padding: '8px 12px',
            borderRadius: '6px',
            zIndex: 10,
            pointerEvents: 'none'
        });
        document.body.appendChild(this.hud);
    }

    _bindEvents() {
        // Mouse drag for pan
        this.canvas.addEventListener('mousedown', e => {
            this.isDragging = true;
            this.lastMouse.x = e.clientX;
            this.lastMouse.y = e.clientY;
            this.simulation.USER_INTERACTION_OCCURING = true;
        });
        window.addEventListener('mousemove', e => {
            if (this.isDragging) {
                const dx = e.clientX - this.lastMouse.x;
                const dy = e.clientY - this.lastMouse.y;
                this.pan.x -= dx * this.zoom;
                this.pan.y -= dy * this.zoom;
                this.lastMouse.x = e.clientX;
                this.lastMouse.y = e.clientY;
                this.updateRenderOptions();
            }
        });
        window.addEventListener('mouseup', () => {
            if (this.isDragging) {
                this.isDragging = false;
                this.simulation.USER_INTERACTION_OCCURING = false;
            }
        });

        // Mouse wheel for zoom (centered on mouse)
        this.canvas.addEventListener('wheel', e => {
            e.preventDefault();
            this.simulation.USER_INTERACTION_OCCURING = true;
            const rect = this.canvas.getBoundingClientRect();
            const mx = e.clientX - rect.left;
            const my = e.clientY - rect.top;
            const worldX = this.pan.x + mx * this.zoom;
            const worldY = this.pan.y + my * this.zoom;
            const zoomFactor = e.deltaY < 0 ? 0.9 : 1.1;
            this.zoom = Math.max(this.minZoom, Math.min(this.maxZoom, this.zoom * zoomFactor));
            this.pan.x = worldX - mx * this.zoom;
            this.pan.y = worldY - my * this.zoom;
            this.updateRenderOptions();

            // Debounce: keep interaction flag true during rapid wheel events
            if (this.wheelDebounceTimer) clearTimeout(this.wheelDebounceTimer);
            this.wheelDebounceTimer = setTimeout(() => {
                this.simulation.USER_INTERACTION_OCCURING = false;
            }, 250);
        }, { passive: false });

        // Optional: handle resize
        window.addEventListener('resize', () => this.updateRenderOptions());
    }

    _bindControlEvents() {
        // Simulation toggle
        const simBtn = document.getElementById('toggle-sim');
        if (simBtn) {
            simBtn.onclick = () => {
                this.simulation.SIM_RUNNING = !this.simulation.SIM_RUNNING;
                simBtn.textContent = this.simulation.SIM_RUNNING ? 'Pause Simulation' : 'Resume Simulation';
            };
        }

        // Symbols toggle
        const symbolsBtn = document.getElementById('toggle-symbols');
        if (symbolsBtn) {
            symbolsBtn.onclick = () => {
                this.simulation.SHOW_SYMBOLS = !this.simulation.SHOW_SYMBOLS;
                symbolsBtn.textContent = this.simulation.SHOW_SYMBOLS ? 'Hide Symbols' : 'Show Symbols';
                // Optionally trigger a render update here
            };
        }

        // Auto zoom toggle
        const autoZoomBtn = document.getElementById('toggle-autozoom');
        if (autoZoomBtn) {
            autoZoomBtn.onclick = () => {
                this.simulation.AUTO_ZOOM = !this.simulation.AUTO_ZOOM;
                autoZoomBtn.textContent = `Auto Zoom: ${this.simulation.AUTO_ZOOM ? 'On' : 'Off'}`;
                if (this.simulation.AUTO_ZOOM && typeof this.autoZoomToParticles === 'function') {
                    this.autoZoomToParticles();
                }
            };
        }
    }

    updateRenderOptions() {
        if (typeof this.simulation.UpdateRenderOptions === 'function') {
            this.simulation.UpdateRenderOptions({
                ulx: this.pan.x,
                uly: this.pan.y,
                inc: this.zoom,
                pxw: this.canvas.width,
                pxh: this.canvas.height
            });
        }
    }

    // Optionally, expose methods to reset or center view
    centerOn(x, y) {
        this.pan.x = x - (this.canvas.width * this.zoom) / 2;
        this.pan.y = y - (this.canvas.height * this.zoom) / 2;
        this.updateRenderOptions();
    }

    setZoom(zoom, centerX, centerY) {
        const worldX = this.pan.x + centerX * this.zoom;
        const worldY = this.pan.y + centerY * this.zoom;
        this.zoom = Math.max(this.minZoom, Math.min(this.maxZoom, zoom));
        this.pan.x = worldX - centerX * this.zoom;
        this.pan.y = worldY - centerY * this.zoom;
        this.updateRenderOptions();
    }

    // Stub for auto zoom (implement as needed)
    autoZoomToParticles() {
        // You can implement logic here to fit all particles in the view.
        // For now, just center on (0,0) and set default zoom:
        this.centerOn(0, 0);
        this.setZoom(1, this.canvas.width / 2, this.canvas.height / 2);
    }

    // Call this from Simulation.Run's callback
    updateHud(fps, realTime, simTime) {
        this.hud.textContent =
            `FPS: ${fps}\n` +
            `Real Time: ${(realTime / 1000).toFixed(2)} s\n` +
            `Sim Time:  ${simTime.toFixed(2)} s`;
    }
}
