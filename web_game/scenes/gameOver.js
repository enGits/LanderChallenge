import { TutorialView } from './tutorial.js';
import { MainMenu } from './menu.js';

export class GameOver {
  constructor(canvas, ctx, switchScene) {
    this.canvas = canvas;
    this.ctx = ctx;
    this.switchScene = switchScene;
    this.stars = Array.from({ length: 100 }, () => ({
      x: Math.random() * this.canvas.width,
      y: Math.random() * this.canvas.height,
    }));
  }

  update() {}

  render() {
    const ctx = this.ctx;
    ctx.fillStyle = 'black';
    ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
    this.stars.forEach(star => {
      ctx.beginPath();
      ctx.arc(star.x, star.y, 1, 0, Math.PI * 2);
      ctx.fillStyle = 'white';
      ctx.fill();
    });

    ctx.textAlign = 'center';
    ctx.fillStyle = 'red';
    ctx.font = '40px Courier New';
    ctx.fillText('Game Over', this.canvas.width / 2, this.canvas.height / 2 - 60);

    ctx.fillStyle = 'white';
    ctx.font = '20px Courier New';
    ctx.fillText('Press ENTER for Main Menu', this.canvas.width / 2, this.canvas.height / 2);
    ctx.fillStyle = 'yellow';
    ctx.fillText('Press T for Tutorial', this.canvas.width / 2, this.canvas.height / 2 + 40);
    // ctx.fillStyle = 'white';
    // ctx.fillText('Press ESC to Quit', this.canvas.width / 2, this.canvas.height / 2 + 80);
  }

  onKeyPress(e) {
    if (e.code === 'Enter' || e.code === 'NumpadEnter') {
      this.switchScene(new MainMenu(this.canvas, this.ctx, this.switchScene));
    } else if (e.code === 'KeyT') {
      this.switchScene(new TutorialView(this.canvas,this.ctx, this.switchScene));
    } 
    // else if (e.code === 'Escape') {
    //   console.log('Exit (placeholder)');
    // }
  }
}
