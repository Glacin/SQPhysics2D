#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include <box2d/box2d.h>

int main() {
    constexpr int layers = 20;
    b2WorldDef world_def = b2DefaultWorldDef();
    world_def.gravity = {0, -10};
    world_def.contactHertz = 64;
    world_def.contactDampingRatio = 8;
    world_def.contactSpeed = 3;
    world_def.restitutionThreshold = 0.1f;
    world_def.enableSleep = false;
    const b2WorldId world = b2CreateWorld(&world_def);

    b2ShapeDef shape_def = b2DefaultShapeDef();
    // Box2D has one Coulomb coefficient; match SQP's static coefficient.
    shape_def.material.friction = 0.75f;
    shape_def.material.restitution = 0;
    shape_def.density = 1;
    b2BodyDef body_def = b2DefaultBodyDef();
    body_def.type = b2_staticBody;
    body_def.position = {0, -0.5f};
    const b2BodyId ground = b2CreateBody(world, &body_def);
    const b2Polygon platform = b2MakeBox(float(layers), 0.5f);
    b2CreatePolygonShape(ground, &shape_def, &platform);

    std::vector<b2BodyId> bodies;
    bodies.reserve(layers * (layers + 1) / 2);
    const b2Polygon box = b2MakeBox(0.5f, 0.5f);
    for (int row = 0; row < layers; ++row) {
        const int count = layers - row;
        const float start_x = -float(count - 1) * 0.5f;
        for (int column = 0; column < count; ++column) {
            body_def = b2DefaultBodyDef();
            body_def.type = b2_dynamicBody;
            body_def.position = {start_x + float(column), 0.5f + float(row)};
            const b2BodyId body = b2CreateBody(world, &body_def);
            b2CreatePolygonShape(body, &shape_def, &box);
            bodies.push_back(body);
        }
    }
    const auto drifts = [&] {
        std::vector<float> result(layers, 0);
        int index = 0;
        for (int row = 0; row < layers; ++row) {
            const int count = layers - row;
            const float start_x = -float(count - 1) * 0.5f;
            for (int column = 0; column < count; ++column) {
                const float x = b2Body_GetPosition(bodies[index++]).x;
                result[row] = std::max(result[row], std::abs(x - start_x - float(column)));
            }
        }
        return result;
    };
    const auto center_sag = [&] {
        std::vector<float> result(layers, 0);
        int index = 0;
        for (int row = 0; row < layers; ++row) {
            const int count = layers - row;
            for (int column = 0; column < count; ++column) {
                const b2Pos position = b2Body_GetPosition(bodies[index++]);
                if (column == (count - 1) / 2 || column == count / 2)
                    result[row] = std::max(result[row], 0.5f + float(row) - position.y);
            }
        }
        return result;
    };
    std::vector<float> half;
    std::vector<float> sag_half;
    for (int tick = 1; tick <= 640; ++tick) {
        b2World_Step(world, 1.0f / 64.0f, 4);
        if (tick == 320) {
            half = drifts();
            sag_half = center_sag();
        }
    }
    const auto final = drifts();
    const auto sag_final = center_sag();
    for (int row = 0; row < layers; ++row)
        std::printf("layer_%d_horizontal_drift_half=%.6f final=%.6f\n", row,
                    double(half[row]), double(final[row]));
    for (int row = 0; row < layers; ++row)
        std::printf("layer_%d_center_sag_half=%.6f final=%.6f\n", row,
                    double(sag_half[row]), double(sag_final[row]));
    b2DestroyWorld(world);
}
