#include "ExtractByPolygon.h"


void PolygonFiller::fillPolygon(QPainter* painter, const QBrush& brush, const int& xsize, const int& ysize, const OGRPolygon *polygon)
{
	QPainterPath outerpath;
	QPainterPath innerpath;
	const OGRLinearRing* exteriorRing = polygon->getExteriorRing();
	QPolygonF outerPolygon;
	for (size_t i = 0; i < exteriorRing->getNumPoints(); i++)
	{
		OGRPoint pt;
		exteriorRing->getPoint(i, &pt);
		outerPolygon.push_back(QPointF(pt.getX(), pt.getY()));

	}
	outerpath.addPolygon(outerPolygon);
	for (size_t n = 0; n < polygon->getNumInteriorRings(); n++)
	{
		QPolygonF innerPolygon;
		const OGRLinearRing* interiorRing = polygon->getInteriorRing(n);
		for (size_t i = 0; i < interiorRing->getNumPoints(); i++)
		{
			OGRPoint pt;
			interiorRing->getPoint(i, &pt);
			innerPolygon.push_back(QPointF(pt.getX(), pt.getY()));
		}
		innerpath.addPolygon(innerPolygon);
		//path.addPolygon(outerPolygon);
	}
	if (polygon->getNumInteriorRings())
		outerpath = outerpath.subtracted(innerpath);
	painter->fillPath(outerpath, brush);

}

RasterMask* PolygonFiller::fillPolygon(const QBrush& br, const float& cellsize, OGRGeometry* poGeometry)
{


	OGRwkbGeometryType tp = poGeometry->getGeometryType();
	OGREnvelope envelope;
	poGeometry->getEnvelope(&envelope);
	int xsize = (int)((envelope.MaxX - envelope.MinX) / cellsize) + 1;
	int ysize = (int)((envelope.MaxY - envelope.MinY) / cellsize) + 1;
	QImage* qImage = new QImage(xsize, ysize, QImage::Format_RGB32);
	qImage->fill(qRgb(0, 0, 0));
	QPainter painter(qImage);
	painter.begin(qImage);
	painter.setRenderHint(QPainter::Antialiasing); //make it look nicer
	painter.setBrush(br);

	QTransform transform = QTransform().translate(-envelope.MinX, -envelope.MaxY) * QTransform().scale((float)xsize / (envelope.MaxX - envelope.MinX), -(float)ysize / (envelope.MaxY - envelope.MinY));
	painter.setTransform(transform);
	if (wkbFlatten(tp) == wkbPolygon)
	{
		fillPolygon(&painter, br, xsize, ysize, (OGRPolygon *)poGeometry);
	}
	else if (wkbFlatten(tp) == wkbMultiPolygon)
	{
		OGRMultiPolygon *polygons = (OGRMultiPolygon*)poGeometry;
		for (size_t n = 0; n < polygons->getNumGeometries(); n++)
		{
			fillPolygon(&painter, br, xsize, ysize, (OGRPolygon*)polygons->getGeometryRef(n));
		}
	}
	painter.end();
	RasterMask* mask = new RasterMask(qImage, envelope.MinX, envelope.MinY, envelope.MaxX, envelope.MaxY);
	return mask;
}

