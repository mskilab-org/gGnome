class Connection extends Base {

  constructor(con) {
    super();
    this.cid = con.cid;
    this.source = con.source ? {sign: Math.sign(con.source), intervalId: Math.abs(con.source)} : null;
    this.sink = con.sink ? {sign: Math.sign(con.sink), intervalId: Math.abs(con.sink)} : null;
    this.title = con.title;
    this.type = con.type;
    this.weight = con.weight;
    this.styleClass = `popovered connection local ${con.type}`;
    this.clipPath = 'url("#clip")';
    this.line = d3.line().curve(d3.curveBasis).x((d) => d[0]).y((d) => d[1]);
  }

  pinpoint(intervalBins) {
    if (this.source) {
      this.source.interval = intervalBins[this.source.intervalId];
      this.source.y = this.source.interval.y;
      this.source.point = this.source.sign > 0 ? this.source.interval.endPoint : this.source.interval.startPoint;
      this.source.place = this.source.sign > 0 ? this.source.interval.endPlace : this.source.interval.startPlace;
      this.touchPlaceX = this.source.place;
      this.touchPlaceY = this.source.y;
      this.touchPlaceSign = this.source.sign;
    }
    if (this.sink) {
      this.sink.interval = intervalBins[this.sink.intervalId];
      this.sink.y = this.sink.interval.y;
      this.sink.point = this.sink.sign > 0 ? this.sink.interval.endPoint : this.sink.interval.startPoint;
      this.sink.place = this.sink.sign > 0 ? this.sink.interval.endPlace : this.sink.interval.startPlace;
      this.touchPlaceX = this.sink.place;
      this.touchPlaceY = this.sink.y;
      this.touchPlaceSign = this.sink.sign;
    }
    this.distance = ((this.source) && (this.sink)) ? d3.format(',')(Math.abs(this.sink.place - this.source.place)) : '-';
  }

  locateInSameFragment(fragment) {
    if (this.source) {
      this.source.scale = fragment.scale;
    }
    if (this.sink) {
      this.sink.scale = fragment.scale;
    }
    this.touchScale = fragment.scale;
    this.identifier = Misc.guid;
  }

  locateInTwoFragments(fragment1, fragment2) {
    this.source.scale = ((this.source.place <= fragment1.domain[1]) && (this.source.place >= fragment1.domain[0])) ? fragment1.scale : fragment2.scale;
    this.sink.scale = ((this.sink.place <= fragment1.domain[1]) && (this.sink.place >= fragment1.domain[0])) ? fragment1.scale : fragment2.scale;
    this.identifier = Misc.guid;
  }

  locateAnchor(fragment) {
    this.kind = 'ANCHOR';
    this.styleClass = `popovered connection anchor`;
    if ((this.source.place <= fragment.domain[1]) && (this.source.place >= fragment.domain[0])) {
      this.source.scale = fragment.scale;
      this.touchPlaceX = this.source.sign > 0 ? this.source.interval.endPlace : this.source.interval.startPlace;
      this.touchPlaceY = this.source.y;
      this.touchPlaceSign = this.source.sign;
      this.fill = d3.rgb(this.sink.interval.color).darker(1);
      this.stroke = '#000';
      this.otherEnd = this.sink;
    } else {
      this.sink.scale = fragment.scale;
      this.touchPlaceX = this.sink.sign > 0 ? this.sink.interval.endPlace : this.sink.interval.startPlace;
      this.touchPlaceY = this.sink.y;
      this.touchPlaceSign = this.sink.sign;
      this.fill = d3.rgb(this.source.interval.color).darker(1);
      this.stroke = '#000';
      this.otherEnd = this.source;
    }
    this.touchScale = fragment.scale;

    this.identifier = Misc.guid;
  }

  get transform() {
    if (this.kind === 'ANCHOR') {
      this.points = [this.touchScale(this.touchPlaceX), this.yScale(this.touchPlaceY)];
      return'translate(' + this.points + ')';
    } else {
      return 'translate(0,0)';
    }
  }

  get render() {
    if (this.kind === 'ANCHOR') {
      this.path = this.arc(this.touchPlaceSign);
    } else {
      this.points = this.type === 'LOOSE' ? this.looseConnectorEndpoints : this.interConnectorEndpoints;
      this.path = this.line(this.points);
    }
    return this.path;
  }

  // Calculate the points for inter-chromosome connections
  get interConnectorEndpoints() {
    var points = [];

    var origin = d3.min([this.source.scale(this.source.place), this.sink.scale(this.sink.place)]);
    var target = d3.max([this.source.scale(this.source.place), this.sink.scale(this.sink.place)]);
    var originSign = (origin === this.source.scale(this.source.place)) ? this.source.sign : this.sink.sign;
    var targetSign = (target === this.source.scale(this.source.place)) ? this.source.sign : this.sink.sign;
    var originY = (origin === this.source.scale(this.source.place)) ? Math.abs(this.source.y) : Math.abs(this.sink.y);
    var targetY = (target === this.source.scale(this.source.place)) ? Math.abs(this.source.y) : Math.abs(this.sink.y);
    var midPointX = 0.5 * origin + 0.5 * target;
    var midPointY = 0.5 * originY + 0.5 * targetY;

    if (this.type === 'ALT') {
      if (Math.abs(this.source.y) === Math.abs(this.sink.y)) {
        points = [
                [origin, this.yScale(originY)],
                [d3.min([origin + Math.sign(originSign) * 5,  midPointX - 5]), this.yScale(originY)],
                [d3.min([origin + Math.sign(originSign) * 25, midPointX - 5]), this.yScale((midPointY + (midPointY < 10 ? 0.5 : 5 )))],
                [midPointX, this.yScale((midPointY + (midPointY < 10 ? 0.75 : 10 )))],
                [d3.max([target + Math.sign(targetSign) * 25, midPointX + 5]), this.yScale((midPointY + (midPointY < 10 ? 0.5 : 5 )))],
                [d3.max([target + Math.sign(targetSign) * 5,  midPointX + 5]), this.yScale(targetY)],
                [target, this.yScale(targetY)]];
      } else {
        points = [
                [origin, this.yScale(originY)],
                [origin + Math.sign(originSign) * 5, this.yScale(originY)],
                [origin + Math.sign(originSign) * 25, this.yScale((originY + Math.sign(targetY - originY) * (originY < 10 ? 0.25 : 5 )))],
                [target + Math.sign(targetSign) * 25, this.yScale((targetY - Math.sign(targetY - originY) * (targetY < 10 ? 0.25 : 5 )))],
                [target + Math.sign(targetSign) * 5, this.yScale(targetY)],
                [target, this.yScale(targetY)]];
      }
    } else {
      points = [
              [origin, this.yScale(originY)],
              [target, this.yScale(targetY)]];
    }
    return points;
  }

  // The array of points forming the loose connections with one endpoint missing
  get looseConnectorEndpoints() {
    return [
      [this.touchScale(this.touchPlaceX), this.yScale(this.touchPlaceY)],
      [this.touchScale(this.touchPlaceX) + this.touchPlaceSign * 15, this.yScale(this.touchPlaceY + (this.touchPlaceY < 10 ? 0.25 : 5 ))],
      [this.touchScale(this.touchPlaceX) + this.touchPlaceSign * 5,  this.yScale(this.touchPlaceY + (this.touchPlaceY < 10 ? 0.75 : 5 ))]];
  }

  // The title for the popover on the connections
  get popoverTitle() {
    return 'Connection #' + this.cid + ' - ' + this.type;
  }

  // The content for the popover on the connections
  get popoverContent() {
    var content = '';
    var array = [
      ['&nbsp;', '<strong>Source</strong>', '<strong>Sink</strong>'], 
      ['Chromosome', ((!this.source) ? 'Unknown' : this.source.interval.chromosome), ((!this.sink) ? 'Unknown' : this.sink.interval.chromosome)], 
      ['Interval', ((!this.source) ? 'Unknown' : (this.source.intervalId + (this.source.sign > 0 ? ' (tail)' : ' (head)'))), ((!this.sink) ? 'Unknown' : (this.sink.intervalId + (this.sink.sign > 0 ? ' (tail)' : ' (head)')))],
      ['Point (chromosome)', ((!this.source) ? 'Unknown' : d3.format(',')(this.source.point)), ((!this.sink) ? 'Unknown' : d3.format(',')(this.sink.point))],
      ['Point (genome)', ((!this.source) ? 'Unknown' : d3.format(',')(this.source.place)), ((!this.sink) ? 'Unknown' : d3.format(',')(this.sink.place))],
      ['Jabba', ((!this.source) ? 'Unknown' : (this.source.y)), ((!this.sink) ? 'Unknown' : (this.sink.y))]
    ];
    array.forEach(function(e,j) {
       content += '<tr><td class="table-label" align="left" width="150" valign="top"><strong>' + e[0] + 
      '</strong></td><td class="table-value" width="100" align="right" valign="top">' + e[1] + 
      '</td><td class="table-value" width="100" align="right" valign="top">' + e[2] + '</td></tr>';
     });
     content += '<tr><td class="table-label" align="left" width="250" valign="top" colspan="2"><strong>Distance</strong></td><td class="table-value" width="100" align="right" valign="top">' + (this.distance) + '</td></tr>';
    return '<div class="row"><div class="col-lg-12"><table width="0" border="0" align="left" cellpadding="0" cellspacing="0"><tbody>' + content + '</tbody></table></div></div>';
  }

  get toString() {
    return `identifier: ${this.identifier},
    cid: ${this.cid},
    source: ${this.source},
    sink: ${this.sink},
    title: ${this.title},
    type: ${this.type}
    weight: ${this.weight}
    `;
  }
}