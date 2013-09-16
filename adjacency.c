//#EID MRR2578
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define buffSize 2048
#define EPSILON .0000001
#define debug 0
struct edge {
  double x1, x2, y1, y2;
  double slope;
  double yIntercept;
};
struct point {
  double xCoord, yCoord;
}; 
struct polygon {
  int number;
  struct point *points;
  int numPoints, maxPoints;
  struct edge *edges;
  int numEdges, maxEdges;
}; 
struct region {
  char *name;
  int number;
  struct polygon *polygons;
  int numPolygons, maxPolygons;
  int knownAdjacent[100];
  int numKnown;
  double minX, minY, maxX, maxY;
}; 
struct map {
  struct region *regions;
  int numRegions, maxRegions;
};
struct map m;

//prototypes
void readFile(FILE *input);
double readPoint(char *line);
void printMap(); 
void printDetails();
void constructEdges();
void compareAll();
int comparePoints(struct point *pointOne, struct point *PointTwo);
int isAdjacentRegions(struct region *regionOne, struct region *regionTwo);
int isAdjacentPolygon(struct polygon *polyOne, struct polygon *polyTwo);
int isIntersecting(struct edge *edgeOne, struct edge *edgeTwo);
double min(double x, double y);
double max(double x, double y);
void freeMap();
//opens file, makes the map, then defines edges, then prints out the map and details of map
int main(int argc, char *argv[]) {
  //open the file
  if(argc > 1) {
    FILE *fp = fopen(argv[1], "r");
    if(fp == NULL)
      printf("No such file exists.\n\n");
    else {
      m.numRegions = 0;
      m.maxRegions = 55;
      m.regions = malloc(55*sizeof(struct region));
      readFile(fp);
      constructEdges();
      compareAll();
      freeMap();
      fclose(fp);
    }
  }
  else
    printf("%s\n", "No file entered.");  
  return 0;
}
//reads in every line, finding regions, polygons, and points
void readFile(FILE *input) {
  char buffer[buffSize];
  char lineDelims[] = {' ', '\t', '\n', ',', '\0'};
  char regionName[buffSize];
  char *line;
  char *previousLine;
  char *result;
  int coordCount = 0;
  int polyCount = 0;
  int regionCount = 0;
  int readingName = 0;
  int regionFlag = 0;
  int polyFlag = 0;
  int pointFlag = 0;
  int lineNum = 0;
  while(line = fgets(buffer, buffSize, input)) {
    if(debug == 1) printf("Line number %d: %s", lineNum, line);
    result = strtok(line, lineDelims);
    while(result != NULL ) {  
      int i;
      for(i = 0; i < (int)strlen(result); i++) {
        if(result[i] == '{') {
          if(debug == 1) printf("Start Region\n");
          if(regionFlag == 1 || polyFlag == 1 || pointFlag == 1) {
            printf("Error at line %d. Improper region start.\n", lineNum);
            exit(0);
          }
          regionFlag = 1;
          readingName = 1;
          if(m.numRegions == m.maxRegions) {
            struct region *temp = realloc(m.regions, (m.maxRegions + 55)*sizeof(struct region));
            m.maxRegions += 55;
            m.regions = temp;
          }  
          if((int)strlen(result) == 1) {
            char *temp = strtok(NULL, lineDelims);
            m.regions[regionCount].name = malloc(((int)strlen(temp) + 1)*sizeof(char));
            strcpy(m.regions[regionCount].name, temp);
          }
          else {
            int c;
            for(c = 1; c < (int)strlen(result); c++)
              regionName[c - 1] = result[c];
            regionName[(int)strlen(result) - 1] = '\0'; 
            m.regions[regionCount].name = malloc(((int)strlen(regionName) + 1)*sizeof(char));  
            strcpy(m.regions[regionCount].name, regionName);
          }        
                  
          m.regions[regionCount].number = regionCount;
          m.regions[regionCount].numKnown = 0;
          m.regions[regionCount].numPolygons = 0;
          m.regions[regionCount].maxPolygons = 5;
          m.regions[regionCount].polygons = malloc(5*sizeof(struct polygon));
          m.regions[regionCount].maxX = -4000000; m.regions[regionCount].minX = 4000000;
          m.regions[regionCount].maxY = -4000000; m.regions[regionCount].minY = 4000000;
          m.numRegions++;
        }
        if(result[i] == '}') {
          if(debug == 1) printf("End Region\n");
          if(polyFlag == 1 || pointFlag == 1 || regionFlag != 1) {
            printf("Error at line %d. Improper region close.\n", lineNum);
            exit(0);
          }
          coordCount = 0;
          polyCount = 0;
          regionFlag = 0;
          regionCount++;
        }
        if(result[i] == '[') {
          if(debug == 1) printf("\tNew Polygon\n");
          if(polyFlag == 1 || pointFlag == 1 || regionFlag != 1) {
            printf("Error at line %d. Improper polygon start.\n", lineNum);
            exit(0);
          }
          readingName = 0;
          polyFlag = 1;
          if(m.regions[regionCount].numPolygons == m.regions[regionCount].maxPolygons) {
            struct polygon *temp = realloc(m.regions[regionCount].polygons, (m.regions[regionCount].maxPolygons + 5)*sizeof(struct polygon));
            m.regions[regionCount].maxPolygons += 5;
            m.regions[regionCount].polygons = temp;
          }
          m.regions[regionCount].polygons[polyCount].number = polyCount;
          m.regions[regionCount].polygons[polyCount].numPoints = 0;
          m.regions[regionCount].polygons[polyCount].maxPoints = 100;
          m.regions[regionCount].polygons[polyCount].points = malloc(100*sizeof(struct point));
          m.regions[regionCount].numPolygons++;
                  
          m.regions[regionCount].polygons[polyCount].numEdges = 0;
          m.regions[regionCount].polygons[polyCount].maxEdges = 100;
          m.regions[regionCount].polygons[polyCount].edges = malloc(100*sizeof(struct edge));
        }
        if(result[i] == ']') {
          if(debug == 1) printf("\tEnd Polygon\n");
          if(polyFlag != 1 || pointFlag == 1 || regionFlag != 1)  {
            printf("Error at line %d. Improper polygon close.\n", lineNum);
            exit(0);
          }
          polyFlag = 0;
          if(comparePoints(&m.regions[regionCount].polygons[polyCount].points[0], &m.regions[regionCount].polygons[polyCount].points[coordCount - 1]) == 0) {
            printf("Error at line %d. Polygon was not closed correctly.\n", lineNum);
            exit(0);
          }
          polyCount++;
          coordCount = 0;
        }
        if(result[i] == '(') {
          if(debug == 1) printf("\t\tX Cord\n");
          if(pointFlag == 1 || polyFlag != 1 || regionFlag != 1)  {
            printf("Error at line %d. Improper x placement.\n", lineNum);
            exit(0);
          }
          pointFlag = 1;
          if(m.regions[regionCount].polygons[polyCount].numPoints == m.regions[regionCount].polygons[polyCount].maxPoints) {
            struct point *temp = realloc(m.regions[regionCount].polygons[polyCount].points, (m.regions[regionCount].polygons[polyCount].maxPoints + 100)*sizeof(struct point));
            m.regions[regionCount].polygons[polyCount].maxPoints += 100;
            m.regions[regionCount].polygons[polyCount].points = temp;
          }
          if(strlen(result) == 1)
            result = strtok(NULL, lineDelims);
          m.regions[regionCount].polygons[polyCount].points[coordCount].xCoord = readPoint(result);
          if(m.regions[regionCount].polygons[polyCount].points[coordCount].xCoord > m.regions[regionCount].maxX)
            m.regions[regionCount].maxX = m.regions[regionCount].polygons[polyCount].points[coordCount].xCoord;

          if(m.regions[regionCount].polygons[polyCount].points[coordCount].xCoord < m.regions[regionCount].minX)
            m.regions[regionCount].minX = m.regions[regionCount].polygons[polyCount].points[coordCount].xCoord;          
        }  
        if(result[i] == ')') {
          if(debug == 1) printf("\t\tY Cord\n");    
          if(pointFlag != 1 || polyFlag != 1 || regionFlag != 1)  {
            printf("Error at line %d. Improper y placement.\n", lineNum);
            exit(0);
          }
          if(strlen(result) == 1)
            m.regions[regionCount].polygons[polyCount].points[coordCount].yCoord = readPoint(previousLine);
          else
            m.regions[regionCount].polygons[polyCount].points[coordCount].yCoord = readPoint(result);

          if(m.regions[regionCount].polygons[polyCount].points[coordCount].yCoord > m.regions[regionCount].maxY)
            m.regions[regionCount].maxY = m.regions[regionCount].polygons[polyCount].points[coordCount].yCoord;  

          if(m.regions[regionCount].polygons[polyCount].points[coordCount].yCoord < m.regions[regionCount].minY)
            m.regions[regionCount].minY = m.regions[regionCount].polygons[polyCount].points[coordCount].yCoord;            

          m.regions[regionCount].polygons[polyCount].numPoints++;    
          coordCount++;
          pointFlag = 0;
        }
      }
      //for regions with spaces in names
      if(readingName == 1 && result[0] != '{') {
        char temp[(int)strlen(result) + (int)strlen(m.regions[regionCount].name) + 2];
        int x;
        for(x = 0; x < (int)strlen(m.regions[regionCount].name); x++)
          temp[x] = m.regions[regionCount].name[x];
        temp[(int)strlen(m.regions[regionCount].name)] = ' ';
        for(x = (int)strlen(m.regions[regionCount].name) + 1; x < (int)strlen(result) + (int)strlen(m.regions[regionCount].name) + 1; x++)
          temp[x] = result[x - (int)strlen(m.regions[regionCount].name) - 1];
        temp[(int)strlen(result) + (int)strlen(m.regions[regionCount].name) + 1] = '\0';
        m.regions[regionCount].name = malloc(((int)strlen(temp) + 1)*sizeof(char));
        strcpy(m.regions[regionCount].name, temp);
      }          
      previousLine = result;
      result = strtok(NULL, lineDelims);
    }
    lineNum++;
  }  
}
//compares two points, ensures closing polygon
int comparePoints(struct point *pointOne, struct point *pointTwo) {
  printf("%1.13f == %1.13f = %d\n", pointOne->xCoord, pointTwo->xCoord, pointOne->xCoord == pointTwo->xCoord);
  printf("%1.13f == %1.13f = %d\n\n", pointOne->yCoord, pointTwo->yCoord, pointOne->yCoord == pointTwo->yCoord);
  if(pointOne->xCoord != pointTwo->xCoord)
    return 0;
  if(pointOne->yCoord != pointTwo->yCoord)
    return 0;
  return 1;
}
//takes a string, either with just a number or (<number> or <number>) and extracts the number
double readPoint(char *line) {
    int i;
  double x;
  int offset = 0;
  char coordString[(int)strlen(line) + 1];
  for(i = 0; i <= (int)strlen(line); i++)
    if(line[i] != ')' && line[i] != '(')
      coordString[i - offset] = line[i];
    else
      offset++;
  x = strtod(coordString, NULL);
  return x;
}
//creates edges for each polygon, edge has a slope, 2 x and 2 y coords, and a y intercept
void constructEdges() {
  int r, p, c;
  int edgeCounter = 0;
  for(r = 0; r < m.numRegions; r++) {
    for(p = 0; p < m.regions[r].numPolygons; p++) {
      edgeCounter = 0;
      for(c = 0; c < m.regions[r].polygons[p].numPoints - 1; c++) {
        if(m.regions[r].polygons[p].numEdges == m.regions[r].polygons[p].maxEdges) {
          struct edge *temp = realloc(m.regions[r].polygons[p].edges, (m.regions[r].polygons[p].maxEdges + 100)*sizeof(struct edge));
          m.regions[r].polygons[p].maxEdges += 100;
          m.regions[r].polygons[p].edges = temp;
        }    
        m.regions[r].polygons[p].edges[edgeCounter].x1 = m.regions[r].polygons[p].points[c].xCoord;
        m.regions[r].polygons[p].edges[edgeCounter].y1 = m.regions[r].polygons[p].points[c].yCoord;      
        m.regions[r].polygons[p].edges[edgeCounter].x2 = m.regions[r].polygons[p].points[c + 1].xCoord;
        m.regions[r].polygons[p].edges[edgeCounter].y2 = m.regions[r].polygons[p].points[c + 1].yCoord;
        if(m.regions[r].polygons[p].edges[edgeCounter].x1 == m.regions[r].polygons[p].edges[edgeCounter].x2)
          m.regions[r].polygons[p].edges[edgeCounter].slope = INFINITY;
        else
          m.regions[r].polygons[p].edges[edgeCounter].slope = (m.regions[r].polygons[p].edges[edgeCounter].y1 - m.regions[r].polygons[p].edges[edgeCounter].y2) / (m.regions[r].polygons[p].edges[edgeCounter].x1 - m.regions[r].polygons[p].edges[edgeCounter].x2);
        
        if(m.regions[r].polygons[p].edges[edgeCounter].slope == 0)
          m.regions[r].polygons[p].edges[edgeCounter].yIntercept = m.regions[r].polygons[p].edges[edgeCounter].y1;
        else if(m.regions[r].polygons[p].edges[edgeCounter].slope == INFINITY)
          m.regions[r].polygons[p].edges[edgeCounter].yIntercept = INFINITY;
        else
          m.regions[r].polygons[p].edges[edgeCounter].yIntercept = m.regions[r].polygons[p].edges[edgeCounter].y1 - (m.regions[r].polygons[p].edges[edgeCounter].slope * m.regions[r].polygons[p].edges[edgeCounter].x1);
        edgeCounter++;        
        m.regions[r].polygons[p].numEdges++;
      }
    }
  }
}
//comparison functions
//sees if two floating points are equal within a certain tolerance
int equalDecimals(double a, double b) {
  if(a == b)
    return 1;
  double n;
  /* find the absolute value of the difference */
    if (a > b)
        n = a-b;
    else 
        n = b-a;
    /* find the larger of the two values */
    double d;
    a = fabs(a);
    b = fabs(b);
    if (a > b)
        d = a;
    else 
        d = b;
    /* now check the ratio of the difference to the larger value */
    if (n/d < EPSILON)
    return 0;
  return 1;
}
//compares all regions to all others to test for adjacency
void compareAll() {
  int r1, r2;
  int firstAdj = 1;
  for(r1 = 0; r1 < m.numRegions; r1++) {
    printf("%s: ", m.regions[r1].name);
    for(r2 = 0; r2 < m.numRegions; r2++)
      if(r1 != r2) {
        if(m.regions[r2].maxX < m.regions[r1].minX)
          continue;
        if(m.regions[r2].minX > m.regions[r1].maxX)
          continue;
        if(m.regions[r2].maxY < m.regions[r1].minY)
          continue;
        if(m.regions[r2].minY > m.regions[r1].maxY)
          continue;
        if(isAdjacentRegions(&m.regions[r1], &m.regions[r2]) == 1)
          if(firstAdj == 1) {
            printf("%s", m.regions[r2].name);
            firstAdj = 0;
          }
          else
            printf(", %s", m.regions[r2].name);        
      }
    printf("\n");
    firstAdj = 1;
  }      
}
//compares a region's polygons to all polygons of other region to test for adjacency
int isAdjacentRegions(struct region *regionOne, struct region *regionTwo) {
  int onePolyCount, twoPolyCount, x;
  //check to see if the region you are comparing is already been found to be adjacent
  for(x = 0; x < regionOne->numKnown; x++)
    if(regionOne->knownAdjacent[x] == regionTwo->number) {
      return 1;
    }  
  for(onePolyCount = 0; onePolyCount < regionOne->numPolygons; onePolyCount++)
    for(twoPolyCount = 0; twoPolyCount < regionTwo->numPolygons; twoPolyCount++)
      if(isAdjacentPolygon(&regionOne->polygons[onePolyCount], &regionTwo->polygons[twoPolyCount]) == 1) {
        regionTwo->knownAdjacent[regionTwo->numKnown] = regionOne->number;
        regionTwo->numKnown++;
        return 1;
      }
  return 0;
}
//compares a polygon's edges to all edges of other polygon to test for adjacency
int isAdjacentPolygon(struct polygon *polyOne, struct polygon *polyTwo) {
  int oneEdgeCount, twoEdgeCount;
  for(oneEdgeCount = 0; oneEdgeCount < polyOne->numEdges; oneEdgeCount++)
    for(twoEdgeCount = 0; twoEdgeCount < polyTwo->numEdges; twoEdgeCount++)
      if(isIntersecting(&polyOne->edges[oneEdgeCount], &polyTwo->edges[twoEdgeCount]) == 1)
        return 1;
  return 0;
}

double min(double x, double y) {
  if(x < y) return x;
  else return y;
}

double max(double x, double y) {
  if(x > y) return x;
  else return y;
}

//tests for adjacency between two edges
int isIntersecting(struct edge *edgeOne, struct edge *edgeTwo) {
  double oneSlope = edgeOne->slope; double twoSlope = edgeTwo->slope;
  double oneX1 = min(edgeOne->x1, edgeOne->x2); double oneY1 = min(edgeOne->y1, edgeOne->y2);
  double oneX2 = max(edgeOne->x1, edgeOne->x2); double oneY2 = max(edgeOne->y1, edgeOne->y2);
  double twoX1 = min(edgeTwo->x1, edgeTwo->x2); double twoY1 = min(edgeTwo->y1, edgeTwo->y2);
  double twoX2 = max(edgeTwo->x1, edgeTwo->x2); double twoY2 = max(edgeTwo->y1, edgeTwo->y2);
  double oneY = edgeOne->yIntercept; double twoY = edgeTwo->yIntercept;
  //two vertical lines
  if(oneSlope == INFINITY && twoSlope == INFINITY) {
    if(oneX1 == twoX1 && oneX1 == twoX2) {
      if(oneY1 == twoY1 && oneY2 == twoY2) return 1;
      if(oneY1 >= twoY1 && oneY2 > twoY2 && equalDecimals(oneY1,twoY2) == 0) return 1;
      if(oneY1 < twoY1 && oneY2 <= twoY2 && equalDecimals(oneY2,twoY1) == 0) return 1;
      if(oneY1 > twoY1 && oneY2 < twoY2) return 1;
      if(oneY1 < twoY2 && oneY2 > twoY1) return 1;
      return 0;
    }
    else return 0;
  }      
  //two horizontal lines
  if(oneSlope == 0.0 && twoSlope == 0.0) {
    if(oneY1 == twoY1 && oneY2 == twoY2) {
      if(oneX1 == twoX1 && oneX2 == twoX2) return 1;
      if(oneX1 >= twoX1 && oneX2 > twoX2 && equalDecimals(oneX1,twoX2) == 0) return 1;
      if(oneX1 < twoX1 && oneX2 <= twoX2 && equalDecimals(oneX2,twoX1) == 0) return 1;
      if(oneX1 > twoX1 && oneX2 < twoX2) return 1;
      if(oneX1 < twoX2 && oneX2 > twoX1) return 1;
      return 0;
    }
    else return 0;
  }
  //two lines with same slope, not horizontal or vertical
  if(equalDecimals(oneSlope, twoSlope) == 1 && oneY == twoY) {
    if(oneX1 == twoX1 && oneX2 == twoX2)
      if(oneY1 == twoY1 && oneY2 == twoY2)
        return 1;
      else return 0;

    if(oneX1 >= twoX1 && oneX2 > twoX2 && equalDecimals(oneX1,twoX2) == 0)
      if(oneY1 >= twoY1 && oneY2 > twoY2 && equalDecimals(oneY1,twoY2) == 0)
        return 1;
      else return 0;

    if(oneX1 < twoX1 && oneX2 <= twoX2 && equalDecimals(oneX2,twoX1) == 0)
      if(oneY1 < twoY1 && oneY2 <= twoY2 && equalDecimals(oneX2,twoX1) == 0)
        return 1;
      else return 0;

    if(oneX1 > twoX1 && oneX2 < twoX2)
      if(oneY1 > twoY1 && oneY2 < twoY2)
        return 1;
      else return 0;

    if(oneX1 < twoX2 && oneX2 > twoX1)
      if(oneY1 < twoY2 && oneY2 > twoY1)
        return 1;
      else return 0;
    return 0;
  }
  return 0;
} 
void freeMap() {
  int r, p;
  for(r = 0; r < m.numRegions; r++) {
    for(p = 0; p < m.regions[r].numPolygons; p++) {
      free(m.regions[r].polygons[p].points);
      free(m.regions[r].polygons[p].edges);
    }
    free(m.regions[r].polygons);
  }
  free(m.regions);
}
//debug functions
//prints the entire map
void printMap() {
  int r, p, c, e;
  for(r = 0; r < m.numRegions; r++) {
    printf("Region #%d: %s\n", r, m.regions[r].name);
    for(p = 0; p < m.regions[r].numPolygons; p++) {
      printf("\tPolygon #%d\n", p);
      for(c = 0; c < m.regions[r].polygons[p].numPoints; c++) {
        printf("\t\tCoordinate #%d: %1.13f, %1.13f\n", c, m.regions[r].polygons[p].points[c].xCoord, m.regions[r].polygons[p].points[c].yCoord);
      }
      for(e = 0; e < m.regions[r].polygons[p].numEdges; e++) {
        printf("\t\tEdge %d: \n\t\t\t%f = x1, %f = y1\n\t\t\t%f = x2, %f = y2\n\t\t\tslope = %f\n\t\t\ty intercept = %f\n", e, m.regions[r].polygons[p].edges[e].x1, m.regions[r].polygons[p].edges[e].y1, m.regions[r].polygons[p].edges[e].x2, m.regions[r].polygons[p].edges[e].y2, m.regions[r].polygons[p].edges[e].slope, m.regions[r].polygons[p].edges[e].yIntercept);
      }
    }
  }
}
//prints simplified version of map
void printDetails() {
  int r, p;
  printf("Num Regions: %d\n", m.numRegions);
  printf("Max Regions: %d\n\n", m.maxRegions);
  for(r = 0; r < m.numRegions; r++) {
    printf("Region #%d: %s\n", r, m.regions[r].name);
    printf("Num Polygons: %d\n", m.regions[r].numPolygons);
    printf("Max Polygons: %d\n", m.regions[r].maxPolygons);
    printf("Max X: %1.13f\n", m.regions[r].maxX);
    printf("Max Y: %1.13f\n", m.regions[r].maxY);
    printf("Min X: %1.13f\n", m.regions[r].minX);
    printf("Min Y: %1.13f\n", m.regions[r].minY);
    for(p = 0; p < m.regions[r].numPolygons; p++) {
      printf("\tPolygon #%d\n", p);
      printf("\t\tNum Points: %d\n", m.regions[r].polygons[p].numPoints);
      printf("\t\tMax Points: %d\n", m.regions[r].polygons[p].maxPoints);
      printf("\t\tNum Edges: %d\n", m.regions[r].polygons[p].numEdges);
      printf("\t\tMax Edges: %d\n", m.regions[r].polygons[p].maxEdges);
    }
  }
}
