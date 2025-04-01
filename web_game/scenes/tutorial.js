import { MainMenu } from './menu.js';
import { OrbitGame } from './orbitGame.js';

export class TutorialView {
  constructor(canvas, ctx, switchScene, previousView = null) {
    this.canvas = canvas;
    this.ctx = ctx;
    this.switchScene = switchScene;
    this.previousView = previousView;
    this.lines = [
      'Press ESC to return',
      'Press ENTER to restart',
      '',
      '## Arrow Keys',
      '- LEFT: Rotate left',
      '- RIGHT: Rotate right',
      '- UP: Increase thrust',
      '- DOWN: Decrease thrust',
      '',
      '## Other',
      '- SPACE: Change focus',
      '- A: Autopilot',
      '- F5: Debug mode',
      '',
      'Press ESC to return'
    ];
  }

  update() {}

  render() {
    const ctx = this.ctx;
    ctx.fillStyle = 'black';
    ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.textAlign = 'center';

    ctx.font = '28px Courier New';
    ctx.fillStyle = 'yellow';
    ctx.fillText('Lander Challenge', this.canvas.width / 2, 100);
    ctx.font = '20px Courier New';
    ctx.fillStyle = 'white';
    ctx.fillText('enGits Lander Challenge', this.canvas.width / 2, 130);

    let y = 180;
    this.lines.forEach(line => {
      ctx.fillStyle = line.startsWith('##') ? 'yellow' : 'white';
      ctx.fillText(line.replace('## ', ''), this.canvas.width / 2, y);
      y += 26;
    });
  }

  onKeyPress(e) {
    if (e.code === 'Escape') {
      if (this.previousView) {
        this.previousView.paused = false; // Resume the game
        this.switchScene(this.previousView);
      } else {
        this.switchScene(new MainMenu(this.canvas, this.ctx, this.switchScene));
      }
    } else if (e.code === 'Enter' || e.code === 'NumpadEnter') {
      this.switchScene(new OrbitGame(this.canvas, this.ctx, this.switchScene, 1));
    }
  }
}
