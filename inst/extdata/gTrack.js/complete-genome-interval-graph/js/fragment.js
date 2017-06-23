class Fragment extends Base {

  constructor(brush) {
    super();
    this.id = Misc.guid;
    this.brush = brush;
    this.selection = null;
    this.domain = null;
    this.range = null;
    this.panelWidth = null;
  }

  get toString() {
    return `fragment: ${this.id}, 
      selection: [${this.selection.join(',')}], 
      domain: [${this.domain.join(',')}]
      range: [${this.range && this.range.join(',')}]
      panelWidth: ${this.panelWidth}`;
  }
}