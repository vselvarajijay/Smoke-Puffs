#include "Fluid.h"
#include <math.h>

void CreateObstacles(Surface dest, int width, int height)
{
    glBindFramebuffer(GL_FRAMEBUFFER, dest.FboHandle);
    glViewport(0, 0, width, height);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    GLuint program = CreateProgram("Fluid.Vertex", 0, "Fluid.Fill");
    glUseProgram(program);

    const int DrawBorder = 1;
    if (DrawBorder) {
        #define T 0.9999f
        float positions[] = { -T, -T, T, -T, T,  T, -T,  T, -T, -T };
        #undef T
        GLuint vbo;
        GLsizeiptr size = sizeof(positions);
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, size, positions, GL_STATIC_DRAW);
        GLsizeiptr stride = 2 * sizeof(positions[0]);
        glEnableVertexAttribArray(PositionSlot);
        glVertexAttribPointer(PositionSlot, 2, GL_FLOAT, GL_FALSE, stride, 0);
        glDrawArrays(GL_LINE_STRIP, 0, 5);
        glDeleteBuffers(1, &vbo);
    }

    const int DrawCircle = 1;
    if (DrawCircle) {
        const int slices = 64;
        float positions[slices*2*3];
        float twopi = 8*atan(1.0f);
        float theta = 0;
        float dtheta = twopi / (float) (slices - 1);
        float* pPositions = &positions[0];
        for (int i = 0; i < slices; i++) {
            *pPositions++ = 0;
            *pPositions++ = 0;

            *pPositions++ = 0.25f * cos(theta) * height / width;
            *pPositions++ = 0.25f * sin(theta);
            theta += dtheta;

            *pPositions++ = 0.25f * cos(theta) * height / width;
            *pPositions++ = 0.25f * sin(theta);
        }
        GLuint vbo;
        GLsizeiptr size = sizeof(positions);
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, size, positions, GL_STATIC_DRAW);
        GLsizeiptr stride = 2 * sizeof(positions[0]);
        glEnableVertexAttribArray(PositionSlot);
        glVertexAttribPointer(PositionSlot, 2, GL_FLOAT, GL_FALSE, stride, 0);
        glDrawArrays(GL_TRIANGLES, 0, slices * 3);
        glDeleteBuffers(1, &vbo);
    }

    // Cleanup
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &vao);

}
