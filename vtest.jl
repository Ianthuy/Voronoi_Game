
using Test

include("vorstellung 1.jl") 

@testset "Test der Kernfunktionen des Voronoi-Spiels" begin

    @testset "Erstellung der Datenstruktur" begin
        p1 = Point(0.0, 0.0)
        p2 = Point(1.0, 0.0)
        p3 = Point(0.0, 1.0)
        
        @test p1.x == 0.0 && p1.y == 0.0
        
        e1 = Edge(p1, nothing, nothing, nothing, nothing)
        @test e1.origin == p1
        
        tri = Dreieck(e1)
        @test tri.edge == e1
        
        gp = GamePoint(10, 20, 1)
        @test gp.owner == 1 && gp.point.x == 10.0
    end

    @testset "Kernalgorithmus von Delaunay" begin
        # --- 测试 point_in_triangle ---
        a = Point(0, 0)
        b = Point(4, 0)
        c = Point(2, 4)
        
        e_ab = Edge(a, nothing, nothing, nothing, nothing)
        e_bc = Edge(b, nothing, nothing, nothing, nothing)
        e_ca = Edge(c, nothing, nothing, nothing, nothing)
        
        e_ab.next = e_bc; e_bc.prev = e_ab
        e_bc.next = e_ca; e_ca.prev = e_bc
        e_ca.next = e_ab; e_ab.prev = e_ca
        
        test_tri = Dreieck(e_ab)
        e_ab.face = test_tri; e_bc.face = test_tri; e_ca.face = test_tri

        p_inside = Point(2, 2)
        p_outside = Point(5, 5)
        p_on_edge = Point(2, 0)

        @test point_in_triangle(p_inside, test_tri) == true
        @test point_in_triangle(p_outside, test_tri) == false
        @test point_in_triangle(p_on_edge, test_tri) == true # 边上的点也算在内

        # --- Test insert_point! ---
        reset_game() 
        delaunay = Delaunay()
        insert_point!(Point(50, 50), delaunay)
        @test length(delaunay.triangles) == 3

        insert_point!(Point(60, 60), delaunay)

        @test length(delaunay.triangles) == 5
    end

    @testset "Berechnung der Voronoi-Diagramm" begin
        reset_game()
        delaunay = Delaunay()
        insert_point!(Point(25, 50), delaunay)
        insert_point!(Point(75, 50), delaunay)
        
        voronoi = compute_voronoi(delaunay)
        @test length(voronoi.edges) > 0
    end

    @testset "Spiel Logik" begin
        reset_game()

        add_point(20, 20) 
        @test length(game_points) == 1
        @test game_points[1].owner == 1
        @test current_player == 2
        @test turn_count == 1

        add_point(80, 80) 
        @test length(game_points) == 2
        @test game_points[2].owner == 2
        @test current_player == 1
        @test turn_count == 2


        remove_nearby_point(20.1, 20.1) 
        @test length(game_points) == 1
        @test game_points[1].owner == 2 


        reset_game()

        add_point(10, 50) 
        add_point(90, 50) 
        

        global game_over = true
        winner_is = determine_winner()
        
    
        @test winner_is === 1


        reset_game()
        add_point(25, 50) 
        add_point(35, 50) 
        winner_is_2 = determine_winner()

        @test winner_is_2 == 2 

        # ---  reset_game ---
        reset_game()
        @test isempty(game_points)
        @test current_player == 1
        @test turn_count == 0
        @test game_over == false
        @test winner === nothing
        @test delaunay_triangulation.bounding_triangle === nothing
    end

end

println("\n Alle Tests wurden abgeschlossen.")

