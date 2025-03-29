import { GameOver } from './gameOver.js';
import { TutorialView } from './tutorial.js';
// import { Body } from '../core/body.js';
// import { Spacecraft } from '../core/spacecraft.js';
// import { Mountain } from '../core/mountain.js';

export class OrbitGame {
  constructor(canvas, ctx, switchScene, level = 1) {
    this.canvas = canvas;
    this.ctx = ctx;
    this.switchScene = switchScene;

    // Game state
    this.level = level;
    this.paused = false;
    this.gameRunning = true;
    this.simTime = 0;
    this.sim_dt = 0.1;
    this.timeFactor = 1.0;

    // Star background
    this.stars = Array.from({ length: 100 }, () => ({
      x: Math.random() * canvas.width,
      y: Math.random() * canvas.height,
      alpha: Math.random(),
      speed: 0.005 + Math.random() * 0.01
    }));

    // Game entities (placeholders)
    this.bodyList = [];   // You'll need to port Body class
    this.spriteList = []; // Manual sprite management for now

    this.showCrashed = false;

    // Placeholder Moon body
    this.planet = {
      radius: 1737400,
      mass: 7.34767309e22,
      name: 'Moon',
      scaleFactor: 0.25 * Math.min(canvas.width, canvas.height) / 1737400,
      color: 'gray'
    };
    this.reference = this.planet;

    // TODO: instantiate Spacecraft, Lander, etc.
    // this.spacecraft = new Spacecraft("craft_01.png", 0.1, this);
    // this.spacecraft.y = this.planet.radius + 300000;
    // this.spacecraft.u = Math.sqrt(Body.G * this.planet.mass / this.spacecraft.y);
    // this.spacecraft.scaleFactor = 10 * this.planet.scaleFactor;

    // this.bodyList.push(this.planet, this.spacecraft);

    // Level setup
    if (this.level === 2) this.setLevel2();
    if (this.level === 3) this.setLevel3();
  }

  setLevel2() {
    // Placeholder: define lander and modify control craft
  }

  setLevel3() {
    // Placeholder
  }

  update(deltaTime) {
    if (this.paused || !this.gameRunning) return;

    // Twinkle stars
    this.stars.forEach(star => {
      star.alpha += star.speed;
      if (star.alpha > 1 || star.alpha < 0) {
        star.speed *= -1;
        star.alpha = Math.max(0, Math.min(1, star.alpha));
      }
    });

    // Update physics simulation
    this.simTime += deltaTime * this.timeFactor;

    // TODO: run physics steps for all bodies
    // this.bodyList.forEach(body => body.updatePhysics(dt));
  }

  render() {
    const ctx = this.ctx;
    const canvas = this.canvas;

    ctx.fillStyle = 'black';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Draw twinkling stars
    this.stars.forEach(star => {
      ctx.beginPath();
      ctx.globalAlpha = star.alpha;
      ctx.fillStyle = 'white';
      ctx.arc(star.x, star.y, 1.5, 0, Math.PI * 2);
      ctx.fill();
    });
    ctx.globalAlpha = 1;

    // HUD (top-left info)
    ctx.save();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.textAlign = 'left';
    ctx.font = '14px Courier New';
    ctx.fillStyle = 'gray';
    let xOffset = canvas.width - 200;
    let yOffset = 30;
    const lineHeight = 20;
    ctx.fillText(`Scale: ${this.reference.scaleFactor.toExponential(2)}`, xOffset, yOffset);
    yOffset += lineHeight;
    ctx.fillText(`Time Factor: ${this.timeFactor.toFixed(2)}`, xOffset, yOffset);
    yOffset += lineHeight;
    ctx.fillText(`Time: ${this.simTime.toFixed(2)}s`, xOffset, yOffset);
    yOffset += lineHeight;
    ctx.fillText(`Reference: ${this.reference.name}`, xOffset, yOffset);
    yOffset += lineHeight;
    ctx.fillText(`Level: ${this.level}`, xOffset, yOffset);
    ctx.restore();

    // Instructions
    ctx.fillStyle = 'white';
    ctx.font = '16px Courier New';
    ctx.fillText('Press T for Tutorial, ESC to Quit', xOffset, canvas.height - yOffset);

    if (this.showCrashed) {
      ctx.font = '40px Courier New';
      ctx.fillStyle = 'red';
      ctx.fillText('CRASHED!', canvas.width / 2, canvas.height / 2);
    }
  }

  onKeyPress(e) {
    switch (e.code) {
      case 'KeyT':
        this.paused = true;
        this.switchScene(new TutorialView(this.canvas, this.ctx, this.switchScene, this));
        break;
      case 'Escape':
        this.gameRunning = false;
        this.triggerGameOver();
        break;
      case 'NumpadAdd':
        this.timeFactor *= 1.1;
        break;
      case 'NumpadSubtract':
        this.timeFactor = Math.max(0.1, this.timeFactor / 1.1);
        break;
    }
  }

  triggerGameOver() {
    this.showCrashed = true;
    setTimeout(() => {
      this.switchScene(new GameOver(this.canvas, this.ctx, this.switchScene));
    }, 4000);
  }
}
